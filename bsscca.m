function [W,r] = bsscca(X,delay)
% bsscca() - Blind Source Separation through Canonical Correlation Analysis
%
% Usage:
%   >> [W,r] = bsscca(X,delay)
%
% Inputs:
%   X     - data matrix (dxN, data channels are rowwise)
%   delay - delay at which the autocorrelation of the sources will be
%           maximized (def: 1)
%
% Output:
%   W     - separation matrix
%   r     - autocorrelation of the estimated sources at the given delay
%
% See also:
%   BSSCCA_IFC, AUTOBSS

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

if nargin < 2, delay = 1; end
if nargin < 1, 
    help cca;
    return;
end


[d,T] = size(X);

% correlation matrices
Y = X(:,delay+1:end);
X = X(:,1:end-delay);
Cyy = (1/T)*Y*Y';
Cxx = (1/T)*X*X';
Cxy = (1/T)*X*Y';
Cyx = Cxy';
invCyy = pinv(Cyy);

% calculate W
[W,r] = eig(pinv(Cxx)*Cxy*invCyy*Cyx);
r = sqrt(abs(real(r)));
[r,I] = sort(diag(r),'descend');
W = W(:,I)';