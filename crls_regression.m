function [Y,H,Hh] = crls_regression(X, opt)
% crls_regression() - Performs automatic EOG artifact correction using
% multiple adaptive regression. The adaptation is made using the Conventional
% Recursive Least Squares Algorithm (CRLS) [1,2]. A forgetting factor can be
% used for dealing with time-varying scenarios. A stable (but slower) version
% of this algorithm is implemented in scrls_regression().
%
% Usage:
%   >>  [Y,H,Hh] = crls_regression( X, opt)
%
% Inputs:
%   X               - Input data matrix
%   opt             - Analysis options (see below)
%   opt.refdata     - Reference signal(s) (dref x N) (def: [])
%   opt.M           - Order of the adaptive filter (def: 3)
%   opt.lambda      - Forgetting factor (def: 0.9999)
%   opt.sigma       - Parameter described in [1] (def: 0.01)
% 
% Outputs:
%   Y   - Output data matrix (artifact corrected)
%   H   - Final filter weights (M*dref x d)
%   Hh  - filter weights evolution (M*dref x d x N)
%
% Notes:
%   1) This function implements the algorithm described in [1].
%
% References:
% [1] P. He et al., Med. Biol. Comput. 42 (2004), 407-412
% [2] S. Haykin. Adaptive Filter Theory, (1996), Prentice Hall
%
% Author: German Gomez-Herrero
%         german.gomezherrero@ieee.org
%         http://www.cs.tut.fi/~gomezher/index.htm
%         Institute of Signal Processing
%         Tampere University of Technology, 2007
%
% See also:
%   POP_CRLS_REGRESSION, POP_SCRLS_REGRESSION, SCRLS_REGRESSION, EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

% overflow value
% ----------------------------
OVERFLOW = 1e12;

if nargin < 1,
    help crls_regression;
    return;
end

if ~exist('opt','var'),
    opt = def_crls_regression;
else
    opt = def_crls_regression(opt);
end

if isempty(opt.refdata),
    error('(crls_regression) I need a reference signal!');
end

sigma = opt.sigma;
lambda = opt.lambda;
M = opt.M;
Xref    = opt.refdata;

[deeg,Leeg] = size(X);
[dref,Lref] = size(Xref);
if Leeg~=Lref, 
    error('(crls_regression) Input and reference data must have the same length'); 
end

% initialization of the adaptation loop
% --------------------------------------
H = zeros(dref*M,deeg);
invR = eye(dref*M)./sigma;
Y = zeros(deeg,Leeg);
if nargout > 2, 
    Hh = zeros(dref*M,deeg,Leeg);
    Hh(:,:,1:M-1) = repmat(H,[1,1,M-1]); 
end

% adaptation loop
% --------------------------------------
if opt.verbose, fprintf('\n(crls_regression) '); end
for i = M:Leeg
    r = Xref(:,i:-1:(i-M+1));
    r = reshape(r', M*dref,1);
    K = (invR*r)./(lambda+r'*invR*r);
    e0 = X(:,i)'-r'*H;  
    if nargout > 2, Hh(:,:,i) = H; end
    H = H+K*e0;    
    invR = (1/lambda)*invR-(1/lambda)*K*r'*invR; 
    Xupdate = r'*H;
    Y(:,i) = (X(:,i)'-Xupdate)';    
    if ~isempty(find(abs(Xupdate(:))>OVERFLOW, 1)),
        error('(crls_regression) Algorithm became unstable');
    end
    if opt.verbose && ~mod(i,floor(Leeg/10)), fprintf('.'); end
end
if opt.verbose,fprintf('[OK]\n');end


% sub-function to initialize the default parameters
% -------------------------------------------------
function [opt] = def_crls_regression(opt)

if ~exist('opt','var') || isempty(opt) || ~isfield(opt,'refdata'),
    opt.refdata = [];
end
if ~isfield(opt, 'verbose') || isempty(opt.verbose),
    opt.verbose = 1;
end
if ~isfield(opt,'sigma') || isempty(opt.sigma),
    opt.sigma = 0.01;
end
if ~isfield(opt, 'lambda') || isempty(opt.lambda),
    opt.lambda = 0.999;
end
if ~isfield(opt, 'M') || isempty(opt.M),
    opt.M = 3;
end

