function [A,B,Lambda,Eiter] = iDMRG_GS(A,B,Lambda,H_loc,Nkeep,Nstep,varargin)

% < Input >
% Lambda : [vector] singular value at the center
% A,B : [rank-3 tensor] for even/odd nodes
%
%  ->-   A    ->-*->-diag(Lambda)->-*->-    B   ->-
%   1    ^     2   1              2   1     ^    2  
%        |3                                 |3
%
% H_loc : [tensor] MPO of interaction Hamiltonian
%       a rank-4 tensor acting on site n. The order of legs
%       of Hs{n} is bottom-top-left-right, where the bottom (top) leg
%       contracts to the physical leg of bra (ket) tensor.
%
% < Option >
% 'Krylov', .. : [numeric] The maximum dimension of the Krylov subspace to
%       be considered, within the Lanczos method used for updating MPS
%       tensors. That is, a tridiagonal matrix of size n-by-n (at maximal)
%       is diagonalized in each tensor update, with n set by 'Krylov',.. .
%       (Default: 5)
% 'tol', .. : [numeric] The tolerance for elements on the +-1 diagonals
%       (those next to the main diagonal) within the Lanczos method. If an
%       element is smaller than the tolerance, then the Lanczos
%       tridiagonalization procedure stops and only the part of the
%       tridiagonal matrix constructed so far is diagonalized.
%       (Default: 1e-8)
%
% < Output >
% Lambda, Gamma : [1 x 2 cells each] Cell arrays of Lambda and Gamma
%       tensors, repectively.
% Eiter : [(numel(taus) x 2 x 2 matrix] Eiter(m,n,k) is the measured energy
%       for an odd (k = 1) or even (k = 2) bond after odd (n = 1) or
%       even (n = 2) bonds are updated, at the m-th "outer" iteration.
%

tobj = tic2;

nKrylov = 5;

% parse options
while ~isempty(varargin)
    if numel(varargin) < 2
        error('ERR: Option should be set by a pair of option name and value.');
    end
    switch varargin{1}
        case 'Krylov'
            nKrylov = varargin{2};
            varargin(1:2) = [];
        case 'tol'
            tol = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: Unknown input.');
    end
end


Skeep = 1e-8;
tol = 1e-8;

% % % check the integrity of input


if ~all(size(A,3) == [size(H_loc,1), size(H_loc,2)])
    error('ERR: The first or second leg of H_loc has dimension inconsistent with the physical leg of an MPS tensor.');
elseif ~all(size(B,3) == [size(H_loc,1), size(H_loc,2)])
    error('ERR: The first or second leg of H_loc has dimension inconsistent with the physical leg of an MPS tensor.');
end

% % %

% show message
disptime('iDMRG: ground state search');
disptime(['Nkeep = ',sprintf('%i',Nkeep),', # of steps = ',sprintf('%i',Nstep)]);



% ground-state energy for each iteration
Eiter = zeros(1,Nstep);
L = 1;
L = updateLeft(L,3,A(1,:,:),H_loc(:,:,end,:),4,A(1,:,:));
R = 1;
R = updateLeft(R,3,permute(B(:,1,:),[2 1 3]),permute(H_loc(:,:,:,1),[1 2 4 3]),4,permute(B(:,1,:),[2 1 3]));
Lambda_old = Lambda;


for itS = (1:Nstep)
    % update L,R and prepare lambda, A,B (step 2,3)
    % move orthogonality center

    T = contract(diag(Lambda),2,2,B,3,1);
    [A2, LambdaR,V] = svdTr(T,3,[1 3], Nkeep,Skeep);
    A2 = permute(A2,[1 3 2]);
%     LambdaR = contract(diag(LambdaR),2,2,V,2,1);

    T = contract(A,3,2,diag(Lambda),2,1,[1 3 2]);
    [U,LambdaL,B2] = svdTr(T,3,1,Nkeep,Skeep);
%     LambdaL = contract(U,2,2,diag(LambdaL),2,1);


    L = updateLeft(L,3,A,H_loc,4,A);
    L = updateLeft(L,3,A2,H_loc,4,A2);
    R = updateLeft(R,3,permute(B,[2 1 3]),permute(H_loc,[1 2 4 3]),4,permute(B,[2 1 3]));
    R = updateLeft(R,3,permute(B2,[2 1 3]),permute(H_loc,[1 2 4 3]),4,permute(B2,[2 1 3]));
    

    % use the trial wavefunction given in step 4.
%     Lambda_new = contract(LambdaR,2,2,diag(1./Lambda),2,1);
%     Lambda_new = contract(Lambda_new,2,2,LambdaL,2,1);
    Aold = contract(A2,3,2,diag(LambdaR./Lambda.*LambdaL),2,1);
    Aold = contract(Aold,3,3,B2,3,1,[1 3 2 4]);

    % variational update
    [Anew,Eiter(itS)] = eigs_2site_GS (L,H_loc,H_loc,R,Aold,nKrylov,Skeep);

    % update next tensor
    [A,S,B] =  svdTr(Anew,4,[1 3],Nkeep,Skeep);
    Lambda = S/norm(Skeep);
    A = permute(A,[1 3 2]);
    


    % stopping criteria
    if sum((Lambda_old-Lambda).^2)<Nkeep*100*tol
        disp('fixed point reached in it = ',sprintf('%i',itS))
        break

    end

    Lambda_old = Lambda;



end


toc2(tobj,'-v');
end




% function for variationally searching ground state by Lanczos method (forked from DMRG_GS_2site.m) 
function [Anew,Enew] = eigs_2site_GS (Hleft,Hcen1,Hcen2,Hright,Aold,nKrylov,tol)
% < Description >
%
% Anew = eigs_2site_GS (Hleft,Hcen1,Hcen2,Hright,Aold,nKrylov,tol)
%
% Update an MPS tensor acting on two neighboring sites, by solving the
% effective Hamiltonian via the Lanczos method.
%
% < Input >
% Hleft : [rank-3 tensor] Left part of the effective Hamiltonian. Its legs
%       are ordered as bottom-top-right.
% Hcen1, Hcen2 : [rank-4 tensor] Center, local parts of the effective
%       Hamiltonian. They are MPO tensors for the current two sites.
% Hright : [rank-3 tensor] Right part of the effective Hamiltonian. Its
%       legs are ordered as bottom-top-left.
% Aold : [rank-4 tensor] Product of two rank-3 ket tensors. Its legs are
%       ordered as (left bond)-(right bond)-(left physical)-(right
%       physical).
% nKrylov, tol : Options for the Lacnzos method. See the documentation for
%       the parent function for details.
%
% < Output >
% Anew : [rank-4 tensor] The approximate ground state of the effective
%       Hamiltonian obtained by the Lanczos method, using the Krylov
%       subspace with dimension up to "nKrylov".
% Enew : [numeric] Expectation value of the effective Hamiltonian with
%       respect to "Anew".

% % % % TODO (start) % % % %

% The collection of rank-4 tensors as Krylov basis; the 5th dimension is
% for indexing different rank-4 tensors
As = zeros([size(Aold,1),size(Aold,2),size(Aold,3),size(Aold,4),nKrylov]);

% Define the first Krylov vector as the input "Aold"
As(:,:,:,:,1) = Aold/sqrt(abs(contract(conj(Aold),4,(1:4),Aold,4,(1:4))));
% normalize; insert "abs" to avoid numerical noise

alphas = zeros(nKrylov,1); % main diagonal elements
betas  = zeros(nKrylov-1,1); % +-1 diagonal elements
cnt = 0; % counter for Krylov subspace dimension

for itn = (1:nKrylov)
    % "matrix-vector" multiplication
    Amul = contract(Hleft,3,2,As(:,:,:,:,itn),4,1);
    Amul = contract(Amul,5,[2 4],Hcen1,4,[3 2]);
    Amul = contract(Amul,5,[3 5],Hcen2,4,[2 3]);
    Amul = contract(Amul,5,[2 5],Hright,3,[2 3],[1 4 2 3]);

    alphas(itn) = real(contract(conj(As(:,:,:,:,itn)),4,(1:4),Amul,4,(1:4)));
    % insert "real" to avoid numerical noise

    cnt = cnt+1;

    if itn < nKrylov
        % orthogonalize, to get the next Krylov vector
        for it2 = (1:2) % do twice to further reduce numerical noise
            T = contract(conj(As(:,:,:,:,1:itn)),5,(1:4),Amul,4,(1:4));
            T = contract(As(:,:,:,:,1:itn),5,5,T,2,1);
            Amul = Amul - T;
        end
    
        % norm
        Anorm = sqrt(abs(contract(conj(Amul),4,(1:4),Amul,4,(1:4)))); % insert "abs" to avoid numerical noise
        
        if Anorm < tol % for numerical stability
            break;
        end
    
        As(:,:,:,:,itn+1) = Amul/Anorm; % normalize
        betas(itn) = Anorm;
    end
end

Hkrylov = diag(betas(1:cnt-1),-1);
Hkrylov = Hkrylov + Hkrylov' + diag(alphas(1:cnt));

[V,D] = eig(Hkrylov);
[~,minid] = min(diag(D));
Anew = contract(As(:,:,:,:,1:cnt),5,5,V(:,minid),2,1);

% compute the epectation value of the effective Hamiltonian with respect to "Anew"
Amul = contract(Hleft,3,2,Anew,4,1);
Amul = contract(Amul,5,[2 4],Hcen1,4,[3 2]);
Amul = contract(Amul,5,[3 5],Hcen2,4,[2 3]);
Amul = contract(Amul,5,[2 5],Hright,3,[2 3],[1 4 2 3]);
Enew = real(contract(conj(Anew),4,(1:4),Amul,4,(1:4)));
% % % % TODO (end) % % % %

end