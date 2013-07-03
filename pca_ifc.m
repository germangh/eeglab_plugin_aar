function [W,A,tmp] = pca_ifc(X, opt)
% pca_ifc() - Interface function for PCA algorithm.
%
% Usage:
%   >> [W,A] = pca_ifc(X [,opt])
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
%   1) If the maximum spread of eigenvalues is violated, the most redundant
%   mixtures will be discarded in the estimation process.
%   4) IMPORTANT: note that if the unkonwn mixing matrix is not of
%   full-column rank we will have that size(A,1)>size(W,1).
%
% See also:
%   PCA, AUTOBSS

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

if ~exist('opt','var') || ~isfield(opt, 'nbsources') || isempty(opt.nbsources),
    opt.nbsources = size(X,1);
end
if ~isfield(opt, 'eigratio') || isempty(opt.eigratio),
    opt.eigratio = 1e6;
end


% reduce the dimensionality of the data
[W] = pca(X,opt.nbsources,opt.eigratio);

A = pinv(W);
tmp = [];