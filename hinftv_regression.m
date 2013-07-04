function [Y,H,Hh] = hinftv_regression(X, opt)
% hinftv_regression() - Performs automatic EOG artifact correction using
% multiple adaptive regression. The adaptation takes place using the H
% infinity norm time-varying algorithm [1].
%
% Usage:
%   >>  [Y,H,Hh] = hinftv_regression( X, opt)
%
% Inputs:
%   X               - Input data matrix (dxN)
%   opt             - Analysis options:
%   opt.refdata     - Reference signal (s) (dref x N) (default: [])
%   opt.M           - Order of the adaptive filter (default: 5)
%   opt.eta         - Factor reflecting a priori knowledge of how close the
%                     estimated filter weights at t=0 are to their optimal
%                     value at that time instant (default: 5e-3)                       
%   opt.rho         - Factor reflecting a priori knowledge of how rapidly
%                     the filter coefficients vary with time
%                     (default: 1e-5)
%   opt.eps         - Positive constant described in [1] (default: 1.5)
% 
% Outputs:
%   Y   - Output data matrix (dxN) (artifact corrected)
%   H   - Final filter weights (M*dref x d)
%   Hh  - filter weights evolution (M*dref x d x N)
%
% Notes:
%   1) This function implements the H infinity TV algorithm proposed in [1].
%
% References:
% [1] S. Puthusserypady and T. Ratnarajah, IEEE Signal Processing Letters
% 12, 816-819
%
%
%
% See also:
%   POP_HINFTV_REGRESSION, EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

% OVERFLOW
OVERFLOW = 1E12;

if nargin < 1,
    help hinftv_regression;
    return;
end

if ~exist('opt','var'),
    opt = def_hinftv_regression;
else
    opt = def_hinftv_regression(opt);
end

if isempty(opt.refdata),
    error('(hinftv_regression) I need a reference signal!');
end

rho = opt.rho;
eta = opt.eta;
M   = opt.M;
eps = opt.eps;
Xref = opt.refdata;
[deeg,Leeg] = size(X);
[dref,Lref] = size(Xref);
if Leeg~=Lref, 
    error('(hinftv_regression) Input data and reference signal must have the same length'); 
end

% initialization of the adaptation loop
% -----------------------------------------------
H = zeros(dref*M,deeg); 
Y = zeros(deeg,Leeg);
Upsilon0 = rho.*eye(M*dref); 
Pi0 = eta*eye(M*dref);
P = Pi0;

if nargout > 2, 
    Hh = zeros(dref*M,deeg,Leeg); 
    Hh(:,:,1:M-1) = repmat(H,[1,1,M-1]); 
end

% adaptation loop
% -----------------------------------------------
if opt.verbose, fprintf('\n(hinftv_regression) '); end
for i = M:Leeg
    r = Xref(:,i:-1:(i-M+1));
    r = reshape(r', M*dref,1);
    P_tilde = inv(inv(P)-eps^(-2)*r*r');
    g = (P_tilde*r)./(1+r'*P_tilde*r);    
    e0 = X(:,i)'-r'*H;
    H = H+g*e0;
    if nargout > 2, Hh(:,:,i) = H; end
    P = inv(inv(P)+(1-eps^(-2))*r*r') + Upsilon0;
    update = r'*H;
    if ~isempty(find(abs(update(:))>OVERFLOW, 1)),
        error('(hinftv_regression) Algorithm became unstable');
    end
    Y(:,i) = (X(:,i)'-update)';
    if opt.verbose && ~mod(i,floor(Leeg/10)), fprintf('.'); end
end
if opt.verbose,fprintf('[OK]\n');end
return;

% sub-function to initialize the default values
% ----------------------------------------------
function [opt] = def_hinftv_regression(opt)

if ~exist('opt','var') || isempty(opt) || ~isfield(opt,'EOGindex'),
    opt.EOGindex = [];
end
if ~isfield(opt, 'verbose') || isempty(opt.verbose),
    opt.verbose = 1;
end
if ~isfield(opt,'eta') || isempty(opt.eta),
    opt.eta = 5e-3;
end
if ~isfield(opt,'rho') || isempty(opt.rho),
    opt.rho = 1e-5;
end
if ~isfield(opt, 'eps') || isempty(opt.eps),
    opt.eps = 1.5;
end
if ~isfield(opt, 'M') || isempty(opt.M),
    opt.M = 3;
end

