function [Y,H,Hout] = lms_regression(X, opt)
% lms_regression() - Performs automatic EOG artifact correction using
% multiple adaptive regression. The adaptation is made using the Least Mean
% Squares (LMS) algorithm [1].
%
% Usage:
%   >>  [Y,H,Hh] = lms_regression( X, opt)
%
% Inputs:
%   X               - Input data matrix (dxN)
%   opt             - Analysis options:
%   opt.refdata     - Reference signal (s) (dref x N) (default: [])
%   opt.M           - Order of the adaptive filter (default: 5)
%   opt.mu          - Learning rate (default: 2e-4)
% 
% Outputs:
%   Y   - Output data matrix (artifact corrected)
%   H   - Final filter weights (M*dref x d)
%   Hh  - filter weights evolution (M*dref x d x N)
%
% References:
% [1] S. Haykin. Adaptive Filter Theory, (1996), Prentice Hall
%
%
% Author: German Gomez-Herrero         
%         http://www.cs.tut.fi/~gomezher/index.htm
%         Institute of Signal Processing
%         Tampere University of Technology, 2007
%
% See also:
%   POP_LMS_REGRESSION, EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

% OVERFLOW
OVERFLOW = 1E12;

if nargin < 1,
    help lms_regression;
    return;
end

if ~exist('opt','var'),
    opt = def_lms_regression;
else
    opt = def_lms_regression(opt);
end

if isempty(opt.refdata),
    error('(lms_regression) I need a reference signal!');
end

mu   = opt.mu;
M    = opt.M;
Xref = opt.refdata;
[deeg,Leeg] = size(X);
[dref,Lref] = size(Xref);
if Leeg~=Lref, 
    error('(lms_regression) Input data and reference signal must have the same length'); 
end

% Initialization of the adaptation loop
% ---------------------------------------------
H = zeros(dref*M,deeg); 
Y = zeros(deeg,Leeg);
if nargout > 2, 
    Hout = zeros(dref*M,deeg,Leeg);
    Hout(:,:,1:M-1) = repmat(H,[1,1,M-1]); 
end

% Adaptation loop
% ---------------------------------------------
if opt.verbose, fprintf('\n(lms_regression) '); end
for i = M:Leeg
    r = Xref(:,i:-1:(i-M+1));
    r = reshape(r', M*dref,1); 
    eupdate = r'*H;
    e0 = X(:,i)'-eupdate;    
    if ~isempty(find(abs(eupdate(:))>OVERFLOW, 1)),
        error('(lms_regression) Algorithm became unstable');
    end
    H = H+mu*r*e0;
    if nargout > 2, Hout(:,:,i) = H; end
    Y(:,i) = e0; 
    if opt.verbose && ~mod(i,floor(Leeg/10)), fprintf('.'); end
end
if opt.verbose,fprintf('[OK]\n');end


% sub-function to initialize the default values
% ----------------------------------------------
function [opt] = def_lms_regression(opt)

if ~exist('opt','var') || isempty(opt) || ~isfield(opt,'refdata'),
    opt.refdata = [];
end
if ~isfield(opt, 'verbose') || isempty(opt.verbose),
    opt.verbose = 1;
end
if ~isfield(opt, 'M') || isempty(opt.M),
    opt.M = 3;
end
if ~isfield(opt, 'mu') || isempty(opt.mu),
    opt.mu = 2e-4;
end

