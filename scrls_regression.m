function [Y,theta,Hh] = scrls_regression(X, opt)
% scrls_regression() - Performs automatic EOG artifact correction using
% multiple adaptive regression. The adaptation is made using the Conventional
% Recursive Least Squares Algorithm (CRLS) [1,2]. A forgetting factor can be
% used for dealing with time-varying scenarios. The precision of the
% computations is set so that stability is guaranteed [3]. If stability is
% not an issue crls_regression() does the same job much faster.
%
% Usage:
%   >>  [Y,H,Hh] = scrls_regression( X, opt)
%
% Inputs:
%   X               - Input data matrix (dxN)
%   opt             - Analysis options (see below)
%   opt.refdata     - Reference signal(s) (dref x N) (default: [])
%   opt.M           - Order of the adaptive filter (def: 3)
%   opt.lambda      - Forgetting factor (def: 0.9999)
%   opt.sigma       - Initialization constant sigma<<1 (def: 0.01)
%   opt.prec        - precision (in bits) to use for the computations (def: 50)
%  
%
% Outputs:
%   Y   - Output data matrix (artifact corrected)
%   H   - Final filter weights (M*dref x d)
%   Hh  - filter weights evolution (M*dref x d x N)
%
% References:
% [1] P. He et al., Med. Biol. Comput. 42 (2004), 407-412
% [2] S. Haykin. Adaptive Filter Theory, (1996), Prentice Hall
% [3] A. P. Liavas and P. A. Regalia, IEEE Trans. Sig. Proc. 47 (1999),
% 88-96
%
%
% See also:
%   POP_SCRLS_REGRESSION, POP_CRLS_REGRESSION, CRLS_REGRESSION, EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

% overflow value for the error
% -----------------------------------
OVERFLOW = 1e12;

if nargin < 1,
    help scrls_regression;
    return;
end

if ~exist('opt','var'),
    opt = def_scrls_regression;
else
    opt = def_scrls_regression(opt);
end

if isempty(opt.refdata),
    error('(scrls_regression) I need a reference signal!');
end

sigma   = opt.sigma;
lambda  = opt.lambda;
M       = opt.M;
Xref    = opt.refdata;

[deeg,Leeg] = size(X);
[dref,Lref] = size(Xref);
if Leeg~=Lref, 
    error('(scrls_regression) Input and reference data must have the same length'); 
end

% compute the maximum allowed precision (in bits)
% ----------------------------------------------
if isempty(opt.prec),    
    if opt.verbose, fprintf('\n(scrls_regression) computing precision bound'); end
    Phi = 0;
    Rk = eye(dref*M)*sigma; R = 0;
    Pk = eye(dref*M)./sigma; P = 0;
    K = 0;
    for i = M:L,
        phi = Xref(:,i:-1:(i-M+1));
        phi = reshape(phi', M*dref,1);
        Phi = max(Phi,norm(phi,1));
        Rk = lambda*Rk+phi*phi';
        R = max(norm(Rk,1),R);
        rk = lambda+phi'*Pk*phi;
        Pk = (1/lambda)*(Pk-(Pk*phi*phi'*Pk)/rk);
        P = max(norm(Pk,1),P);
        K = max(norm(Pk,1)*norm(Rk,1),K);
        if opt.verbose && ~mod(i,floor(Leeg/10)), fprintf('.'); end
    end
    A1=Phi^2*(K+P*Phi^2)^2;
    A2=lambda*(1-lambda)*Phi^2;
    epsilon = (1-lambda)/(K^2*P);
    if isinf(A1^2+A1*A2),
        k = kmax;
    else
        p1=(lambda/Phi^2)*(1-sqrt((A1^2+A1*A2)/(A1+A2)^2));
        epsilon = epsilon*(p1-(A1*p1^2)/(lambda*(1-lambda)*(lambda-Phi^2*p1)));
        k = log2(1/epsilon);
    end
    if opt.verbose,fprintf('[OK]\n');end
else
    k = opt.prec;
end

% initialization of the adaptation loop
% -------------------------------------------
theta = zeros(dref*M,deeg);
P = eye(dref*M)./sigma;
Y = zeros(deeg,Leeg);

% filter recursion
% -------------------------------------------
if opt.verbose, fprintf('\n(scrls_regression) '); end
if nargout < 2, 
    Hh = zeros(dref*M,deeg,Leeg);
    Hh(:,:,1:M-1) = repmat(theta,[1,1,M-1]); 
end
for i = M:Leeg
    phi = Xref(:,i:-1:(i-M+1));
    phi = reshape(phi', M*dref,1);
    u = X(:,i);
    e = u'-phi'*theta;
    e = fl(e,k);
    r = lambda + phi'*P*phi;
    r = fl(r,k);
    v = ((P*phi));
    v = fl(v,k);
    v = v/r;
    v = fl(v,k);
    % update filter weights
    theta = theta + repmat(v,1,deeg).*repmat(e,M*dref,1);
    theta = fl(theta,k);
    if nargout > 2, Hh(:,:,i) = theta; end
    % update inverse of correlation matrix
    Pnew = P*phi;
    Pnew = fl(Pnew,k);
    Pnew = Pnew*phi';
    Pnew = fl(Pnew,k);
    Pnew =Pnew*P;
    Pnew = fl(Pnew,k);
    Pnew = Pnew./r;
    Pnew = fl(Pnew,k);
    Pnew = P - Pnew;
    m = (1/lambda);
    m = fl(m,k);
    Pnew = m*Pnew;
    Pnew = fl(Pnew,k);
    P = Pnew;
    fout = phi'*theta;
    if ~isempty(find(abs(fout(:))>OVERFLOW, 1)),
        error('(scrls_regression) Algorithm became unstable');
    end
    Y(:,i) = (u'-fout)';
    if opt.verbose && ~mod(i,floor(Leeg/10)), fprintf('.'); end
end
if opt.verbose,fprintf('[OK]\n');end


% sub-function to initialize the default parameter values
% ---------------------------------------------------------
function [opt] = def_scrls_regression(opt)

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
if ~isfield(opt, 'prec') || isempty(opt.prec),
    opt.prec = 50;
end


