function [W,A,r] = bsscca_ifc(X, opt)
% cca_ifc() - Interface function for CCA BSS algorithm.
%
% Usage:
%   >> [W,A,r] = bsscca_ifc(X [,opt])
% 
% Inputs:
%   X               - data matrix (dxN, data channels are rowwise)
%   opt.nbsources   - number of sources (def: same as number of data channels)
%   opt.delay       - delay to use for the CCA estimation (def: 1)
%   opt.eigratio    - maximum spread of data covariance eigenvalues. The
%                     spread is measured as lambda_max/lambda_min
%                     (def: 1e6)
% 
% Output:
%   W               - Separation matrix
%   A               - Mixing matrix
%   r               - Autocorrelation of the estimated sources
%
% Notes:
%   1) If the maximum spread of eigenvalues is violated, the most redundant
%   mixtures will be discarded in the estimation process.
%   4) IMPORTANT: note that if the unkonwn mixing matrix is not of
%   full-column rank we will have that size(A,1)>size(W,1).
%
% See also:
%   BSSCCA, AUTOBSS


% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com
 
if ~exist('opt','var') || ~isfield(opt, 'nbsources') || isempty(opt.nbsources),
    opt.nbsources = size(X,1);
end
if ~isfield(opt, 'delay') || isempty(opt.delay),
    opt.delay = 1;
end

if ~isfield(opt, 'eigratio') || isempty(opt.eigratio),
    opt.eigratio = 1e6;
end


% reduce the dimensionality of the data
if opt.eigratio < Inf || opt.nbsources < size(X,1),
    [Wpca,X] = pca(X,opt.nbsources,opt.eigratio);
else
    Wpca = eye(opt.nbsources);
end

[W,r] = bsscca(X, opt.delay);
W = W*Wpca;
A = pinv(W);

