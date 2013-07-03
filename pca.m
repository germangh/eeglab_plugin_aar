function [W,Y] = pca(X,nbc,eigratio)
% PCA principal component analysis
%   [W,Y] = pca(X,NBC,EIGRATIO) returns the PCA matrix W and the principal
%   components Y corresponding to the data matrix X (realizations
%   columnwise). The number of components is NBC components unless the
%   ratio between the maximum and minimum covariance eigenvalue is below
%   EIGRATIO. In such a case, the function will return as few components as
%   are necessary to guarantee that such ratio is greater than EIGRATIO.
% 
% See also:
%   AUTOBSS, EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

if nargin < 3, eigratio = 1e8; end
if nargin < 2, nbc = size(X,1); end


C = cov(X');
[V,D] = eig(C);
[val,I] = sort(abs(diag(D)),'descend');

while val(1)/val(nbc)>eigratio,
    nbc = nbc-1;
end
V = V(:,I(1:nbc));
tmp =diag(D);
D = diag(tmp(I(1:nbc)));
W = D^(-.5)*V';
Y = W*X;
