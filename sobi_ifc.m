function [W,A,out3] = sobi_ifc(X, opt)
% sobi_ifc() - Interface function for SOBI algorithm.
%
% Usage:
%   >> [W,A] = sobi_ifc(X [,opt])
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
%   2) See SOBI for references and credit for the SOBI code.
%   3) If the maximum spread of eigenvalues is violated, the most redundant
%   mixtures will be discarded in the estimation process.
%   4) IMPORTANT: note that if the unkonwn mixing matrix is not of
%   full-column rank we will have that size(A,1)>size(W,1).
%
% See also:
%   SOBI, AUTOBSS


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

% reduce the dimensionality of the data
[Wpca,X] = pca(X,size(X,1),opt.eigratio);

A = sobi(X, opt.nbsources, opt.nbcorr);
W = pinv(A);
A = pinv(Wpca)*A;
W = W*Wpca;
out3 = [];