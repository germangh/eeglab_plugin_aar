function [W,A,out3] = iwasobi_ifc(X, opt)
% sobi_ifc() - Interface function for iWASOBI algorithm.
%
% Usage:
%   >> [W,A] = iwasobi_ifc(X [,opt])
% 
% Inputs:
%   X               - data matrix (dxN, data channels are rowwise)
%   opt.nbsources   - number of sources (def: same as number of data channels)
%   opt.nbcorr      - number of correlation matrices (def: min(100,N/3))
%   opt.eigratio    - maximum spread of data covariance eigenvalues. The
%                     spread is measured as lambda_max/lambda_min
%                     (def: 1e6)
% 
% Output:
%   W               - Separation matrix
%   A               - Mixing matrix
%
% Notes:
%   1) You can use this function as a model to build interfaces to other
%   ICA algorithms. Any interface function accepts the same input
%   parameters and returns the separability matrix.
%   2) See iWASOBI for references and credit for the iWASOBI code.
%   3) If the maximum spread of eigenvalues is violated, the most redundant
%   mixtures will be discarded in the estimation process.
%   4) IMPORTANT: note that if the unkonwn mixing matrix is not of
%   full-column rank we will have that size(A,1)>size(W,1).
%
% See also:
%   IWASOBI, AUTOBSS


% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

if ~exist('opt','var') || ~isfield(opt, 'nbsources') || isempty(opt.nbsources),
    opt.nbsources = size(X,1);
end
if ~isfield(opt, 'nbcorr') || isempty(opt.nbcorr),
    opt.nbcorr = min(100,floor(size(X,2)/3));
end
if ~isfield(opt, 'eigratio') || isempty(opt.eigratio),
    opt.eigratio = 1e6;
end
if ~isfield(opt,'ar_order') || isempty(opt.ar_order),
    opt.ar_order=10;
end
if ~isfield(opt,'rmax') || isempty(opt.rmax),
    opt.rmax = 0.99;
end
if ~isfield(opt,'eps0') || isempty(opt.eps0),
    opt.eps0 = 5e-7;
end

% reduce the dimensionality of the data
[Wpca,X] = pca(X,opt.nbsources,opt.eigratio);

W = iwasobi(X, opt.ar_order,opt.rmax,opt.eps0);
W = W*Wpca;
A = pinv(W);
out3 = [];