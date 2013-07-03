% pop_crls_regression() - Automatic EOG correction using Conventional
% Recursive Least Squares (CRLS) regression
%
% Usage:
%   >> OUTEEG = pop_crls_regression(INEEG, EOGchans, M, lambda, sigma, evchans)
%
% Inputs:
%   INEEG    - input EEG dataset
%   EOGchans - indexes of the EOG reference channels
%   M        - Order of the adaptive filter
%   lambda   - Forgetting factor
%   sigma    - Parameter described in [1]
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
%
% See also:
%   CRLS_REGRESSION, POP_SCRLS_REGRESSION, SCRLS_REGRESSION, EEGLAB

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

function [EEG,com] = pop_crls_regression(EEG, EOGindex, M, lambda,sigma)

com = '';

if nargin < 1,
    help pop_lms_regression;
    return;
end
if nargin < 2, EOGindex = []; end
if nargin < 3, M = 3; end
if nargin < 4, lambda = 0.9999; end
if nargin < 5, sigma = 0.01; end
if nargin < 6, evchans = []; end

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
if nargin < 6,    
    uigeom = {[1.5 1] [1.5 1] [1.5 1] [1.5 1] [1.5 1]};
    uilist = { { 'style' 'text'      'string'    'EOG channel indexes:'} ...
        {'style'  'edit' 'string'    num2str(EOGindex)} ...
        {'style'  'text'      'string'    'Filter order (M):'} ...
        {'style'  'edit'      'string'    num2str(M)} ...
        {'style'  'text'      'string'    'Forgetting factor (lambda):'} ...
        {'style'  'edit'      'string'    num2str(lambda)} ...
        {'style'  'text'      'string'    'Sigma:'} ...
        {'style'  'edit'      'string'    num2str(sigma)} ...
        {'style'  'text'      'string'     'Store filter weights for channels:'} ...
        {'style'  'edit'      'string'    num2str(evchans)} ...        
        };
    guititle = 'Correct EOG using CRLS regression -- pop_crls_regression()';
    result = inputgui( uigeom, uilist, 'pophelp(''pop_crls_regression'')', guititle, [], 'normal');

    if isempty(result), return; end;

    % reading params
    % -------------------
    EOGindex = eval(['[' result{1} ']']);
    M = eval(['[' result{2} ']']);
    lambda = eval(['[' result{3} ']']); 
    sigma = eval(['[' result{4} ']']); 
    evchans = eval(['[' result{5} ']']);
end

% build the options structure
% -----------------------------
opt.refdata = reshape(EEG.data(EOGindex,:,:),length(EOGindex),EEG.pnts*EEG.trials);
opt.M = M;
opt.lambda = lambda;
opt.sigma = sigma;
EEGindex = find(~ismember(1:EEG.nbchan,EOGindex));

% run the EOG correction
% ----------------------
tmp_in = reshape(EEG.data,EEG.nbchan,EEG.pnts*EEG.trials);
index1 = intersect(EEGindex,evchans);
index2 = intersect(EEGindex,find(~ismember(1:EEG.nbchan,evchans)));
if ~isempty(index1),
    [tmp,H,Hh] = crls_regression(tmp_in(index1,:), opt);
    EEG.data(index1,:,:) = reshape(tmp,[length(index1),EEG.pnts,EEG.trials]);
    EEG.Hh = Hh;
    EEG.evchans = evchans;
end
if ~isempty(index2),
    tmp = crls_regression(tmp_in(index2,:), opt);
    EEG.data(index2,:,:) = reshape(tmp,[length(index2),EEG.pnts,EEG.trials]);
end

% command history
% -------------------
com = sprintf( '%s = pop_crls_regression( %s, [%s], %s, %s, %s, [%s]);', inputname(1), ...
    inputname(1), num2str(EOGindex), num2str(M), num2str(lambda),num2str(sigma),num2str(evchans));

return;