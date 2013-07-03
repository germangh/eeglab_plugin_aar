% pop_scrls_regression() - Automatic EOG correction using Stable Conventional
% Recursive Least Squares (SCRLS) regression
%
% Usage:
%   >> OUTEEG = pop_scrls_regression(INEEG, EOGchans, M, lambda, sigma,
%   prec, evchans)
%
% Inputs:
%   INEEG    - input EEG dataset
%   EOGchans - indexes of the EOG reference channels
%   M        - Order of the adaptive filter
%   lambda   - Forgetting factor
%   sigma    - Initialization constant (sigma<<1)
%   prec     - precision (in bits) to use for the computations
%   evchans  - Channel indexes for which the evolution of the filter
%              weights will be stored. This can slow down the
%              computations or even make the program crash because of a
%              memory overflow. Notice that the memory required to store
%              the evolution of K channels is: L*K*N*8 bytes where L is the
%              filter order and N is the number of data samples.
%
% Outputs:
%   OUTEEG  - output dataset
%
% References:
% [1] P. He et al., Med. Biol. Comput. 42 (2004), 407-412
% [2] S. Haykin. Adaptive Filter Theory, (1996), Prentice Hall
% [3] A. P. Liavas and P. A. Regalia, IEEE Trans. Sig. Proc. 47 (1999),
% 88-96
%
% Note:
% In order to guarantee stability of the conventional RLS algorithm we
% impose bounds to the precision of the computations. The current
% implementation does this in a very inefficient way which leads to a
% considerable increase in computation time. Therefore it is advisable to
% use first the standard RLS algorithm and use the stable RLS algorithm
% only if the standard RLS becomes unstable.         
%
% See also:
%   SCRLS_REGRESSION, POP_CRLS_REGRESSION, CRLS_REGRESSION, EEGLAB

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

function [EEG,com] = pop_scrls_regression(EEG, EOGindex, M, lambda,sigma,prec)

com = '';

if nargin < 1,
    help pop_lms_regression;
    return;
end
if nargin < 2, EOGindex = []; end
if nargin < 3, M = 3; end
if nargin < 4, lambda = 0.9999; end
if nargin < 5, sigma = 0.01; end
if nargin < 6, prec = 50; end
if nargin < 7, evchans = []; end

if isempty(EEG.data)
    disp('(pop_crls_regression) error: cannot clean an empty dataset'); 
    return;
end;

% try to guess which are the EOG channels
% --------------------------------------------------
if isempty(EOGindex),
    for i = 1:length(EEG.chanlocs),
        labels = EEG.chanlocs(i).labels;
        if ~isempty(strfind(lower(labels),'eog')),
            EOGindex = [EOGindex i];
        end
    end
end

% display input dialog
% ---------------------------
if nargin < 7,    
    uigeom = {[1.5 1] [1.5 1] [1.5 1] [1.5 1] [1.5 1] [1.5 1]};
    uilist = { { 'style' 'text'      'string'    'EOG channel indexes:'} ...
        {'style'  'edit' 'string'    num2str(EOGindex)} ...
        {'style'  'text'      'string'    'Filter order (M):'} ...
        {'style'  'edit'      'string'    num2str(M)} ...
        {'style'  'text'      'string'    'Forgetting factor (lambda):'} ...
        {'style'  'edit'      'string'    num2str(lambda)} ...
        {'style'  'text'      'string'    'Sigma:'} ...
        {'style'  'edit'      'string'    num2str(sigma)} ...
        {'style'  'text'      'string'    'Precision (in bits):'} ...
        {'style'  'edit'      'string'    num2str(prec)} ...
        {'style'  'text'      'string'     'Store filter weights for channels:'} ...
        {'style'  'edit'      'string'    num2str(evchans)} ...
        };
    guititle = 'Correct EOG using SCRLS regression -- pop_scrls_regression()';
    result = inputgui( uigeom, uilist, 'pophelp(''pop_scrls_regression'')', guititle, [], 'normal');

    if isempty(result), return; end;

    % reading params
    % -------------------
    EOGindex = eval(['[' result{1} ']']);
    M = eval(['[' result{2} ']']);
    lambda = eval(['[' result{3} ']']); 
    sigma = eval(['[' result{4} ']']); 
    prec = eval(['[' result{5} ']']); 
    evchans = eval(['[' result{6} ']']);
end

% build the options structure
% -----------------------------
opt.refdata = reshape(EEG.data(EOGindex,:,:),length(EOGindex),EEG.pnts*EEG.trials);
opt.M = M;
opt.lambda = lambda;
opt.sigma = sigma;
opt.prec = prec;
EEGindex = find(~ismember(1:EEG.nbchan,EOGindex));

% run the EOG correction
% ----------------------
tmp_in = reshape(EEG.data,EEG.nbchan,EEG.pnts*EEG.trials);
index1 = intersect(EEGindex,evchans);
index2 = intersect(EEGindex,find(~ismember(1:EEG.nbchan,evchans)));
if ~isempty(index1),
    [tmp,H,Hh] = scrls_regression(tmp_in(index1,:), opt);
    EEG.data(index1,:,:) = reshape(tmp,[length(index1),EEG.pnts,EEG.trials]);
    EEG.Hh = Hh;
    EEG.evchans = evchans;
end
if ~isempty(index2),
    tmp = scrls_regression(tmp_in(index2,:), opt);
    EEG.data(index2,:,:) = reshape(tmp,[length(index2),EEG.pnts,EEG.trials]);
end

% command history
% -------------------
com = sprintf( '%s = pop_scrls_regression( %s, [%s], %s, %s, %s, %s, [%s]);', inputname(1), ...
    inputname(1), num2str(EOGindex), num2str(M), num2str(lambda), ...
    num2str(sigma),num2str(prec),num2str(evchans));

return;