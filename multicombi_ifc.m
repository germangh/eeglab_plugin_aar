function [W,A,out3] = multicombi_ifc(X, opt)
% multicombi_ifc() - Interface function to MULTICOMBI
%
% Usage:
%   >> [W,A] = efica_ifc(X [,opt])
% 
% Inputs:
%   X               - data matrix (dxN, data channels are rowwise)
%   opt.nbsources   - number of sources (def: same as number of data channels)
%   opt.eigratio    - maximum spread of data covariance eigenvalues. The
%                     spread is measured as lambda_max/lambda_min
%                     (def: 1e6)
% 
% Output:
%   W               - Separation matrix
%   A               - Mixing matrix
%
% Notes:
%   1) See comments in install.txt and multicombi.m for references and credit
%      for the MULTICOMBI code.
%   2) If the maximum spread of eigenvalues is violated, the most redundant
%      mixtures will be discarded in the estimation process.
%   3) IMPORTANT: note that if the unkonwn mixing matrix is not of
%      full-column rank we will have that size(A,1)>size(W,1).
%
% See also:
%   MULTICOMBI, POP_AUTOBSSEMG, POP_AUTOBSSEOG, AUTOBSS, EEGLAB


% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

if ~exist('opt','var') || ~isfield(opt, 'nbsources') || isempty(opt.nbsources),
    opt.nbsources = size(X,1);
end
if ~isfield(opt, 'eigratio') || isempty(opt.eigratio),
    opt.eigratio = 1e6;
end

% reduce the dimensionality of the data
[Wpca,X] = pca(X,opt.nbsources,opt.eigratio);

% call MULTICOMBI
[W] = multicombi(X);
A = inv(W);
A = pinv(Wpca)*A;
W = W*Wpca;
out3 = [];