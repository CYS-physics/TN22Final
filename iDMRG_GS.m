function [Lambda,Gamma,Eiter] = iDMRG_GS(Lambda,Gamma,H,Nkeep,Nmax)

% < Input >
% Lambda : [1 x 2 cell] Lambda{1} and Lambda{2} contain the singular values
%       at odd and even bond, respectively, as column vectors. An odd
%       (even) bond sits just on the right of an odd (even) site.
% Gamma : [1 x 2 cell] Gamma{1} and Gamma{2} are rank-3 "Gamma" tensors for
%       odd and even sites, respectively. Their legs are ordered as left-
%       right-physical(bottom).
%       In an infinite MPS, the tensors are repeated as follows (here the
%       numbers next to the legs indicate their orders):
%
% ->-Gamma{1}->-*->-diag(Lambda{1})->-*->-Gamma{2}->-*->-diag(Lambda{2})->- 
% 1    ^     2   1                 2   1     ^    2   1                 2 
%      |3                                    |3
%
% H : [tensor] Two-site interaction Hamiltonian. Its leg convention is as
%       below:
%
%    2         4        [ 2 (4) is to be contracted with the third leg
%    ^         ^          of Gamma{1} (Gamma{2}) ]
%    |   ...   |
%   [     H     ]
%    |   ...   |
%    ^         ^
%    1         3
% < Output >
% Lambda, Gamma : [1 x 2 cells each] Cell arrays of Lambda and Gamma
%       tensors, repectively.
% Eiter : [(numel(taus) x 2 x 2 matrix] Eiter(m,n,k) is the measured energy
%       for an odd (k = 1) or even (k = 2) bond after odd (n = 1) or
%       even (n = 2) bonds are updated, at the m-th "outer" iteration.
%

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