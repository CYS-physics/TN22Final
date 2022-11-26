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

