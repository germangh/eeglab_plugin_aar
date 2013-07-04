% pop_lms_regression() - Automatic EOG correction using Least Mean Squares
% (LMS) regression
%
% Usage:
%   >> [OUTEEG] = pop_lms_regression(INEEG, EOGchans, M, mu, evchans)
%
% Inputs:
%   INEEG    - input EEG dataset
%   EOGchans - indexes of the EOG reference channels
%   M        - Order of the adaptive filter
%   mu       - Learning rate
%   evchans  - Channel indexes for which the evolution of the filter
%              weights will be stored. This can slow down the
%              computations or even make the program crash because of a
%              memory overflow. Notice that the memory required to store
%              the evolution of K channels is: L*K*N*8 bytes where L is the
%              filter order and N is the number of data samples. The
%              evolution of the weights will be stored in the field .Hh of
%              the output EEG structure.
%
% Outputs:
%   OUTEEG  - output dataset
%
% Author: German Gomez-Herrero
%         <german.gomezherrero@ieee.org>
%         http://www.cs.tut.fi/~gomezher/index.htm
%         Institute of Signal Processing
%         Tampere University of Technology, 2007
%
% References:
% [1] S. Haykin. Adaptive Filter Theory, (1996), Prentice Hall
%
% See also:
%   LMS_REGRESSION, EEGLAB

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

function [EEG,com] = pop_lms_regression(EEG, EOGindex, M, mu, evchans)

com = '';

if nargin < 1,
    help pop_lms_regression;
    return;
end
if nargin < 2, EOGindex = []; end
if nargin < 3, M = 3; end
if nargin < 4, mu = 1e-6; end
if nargin < 5, evchans = []; end

if isempty(EEG.data)
    disp('(pop_lms_regression) error: cannot clean an empty dataset'); 
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
if nargin < 5,    
    uigeom = {[1.5 1] [1.5 1] [1.5 1] [1.5 1]};
    uilist = { { 'style' 'text'      'string'    'EOG channel indexes:'} ...
        {'style'  'edit' 'string'    num2str(EOGindex)} ...
        {'style'  'text'      'string'    'Filter order (M):'} ...
        {'style'  'edit'      'string'    num2str(M)} ...
        {'style'  'text'      'string'    'Learning rate (mu):'} ...
        {'style'  'edit'      'string'    num2str(mu)} ...        
        {'style'  'text'      'string'     'Store filter weights for channels:'} ...
        {'style'  'edit'      'string'    num2str(evchans)} ...        
        };
    guititle = 'Correct EOG using LMS regression -- pop_lms_regression()';
    result = inputgui( uigeom, uilist, 'pophelp(''pop_lms_regression'')', guititle, [], 'normal');

    if isempty(result), return; end;

    % reading params
    % -------------------
    EOGindex = eval(['[' result{1} ']']);
    M = eval(['[' result{2} ']']);
    mu = eval(['[' result{3} ']']);    
    evchans = eval(['[' result{4} ']']);
end

% build the options structure
% -----------------------------
opt.refdata = reshape(EEG.data(EOGindex,:,:),length(EOGindex),EEG.pnts*EEG.trials);
opt.M = M;
opt.mu = mu;
EEGindex = find(~ismember(1:EEG.nbchan,EOGindex));

% run the EOG correction
% ----------------------
tmp_in = reshape(EEG.data,EEG.nbchan,EEG.pnts*EEG.trials);
index1 = intersect(EEGindex,evchans);
index2 = intersect(EEGindex,find(~ismember(1:EEG.nbchan,evchans)));
if ~isempty(index1),
    [tmp,H,Hh] = lms_regression(tmp_in(index1,:), opt);
    EEG.data(index1,:,:) = reshape(tmp,[length(index1),EEG.pnts,EEG.trials]);
    EEG.Hh = Hh;
    EEG.evchans = evchans;
end
if ~isempty(index2),
    tmp = lms_regression(tmp_in(index2,:), opt);
    EEG.data(index2,:,:) = reshape(tmp,[length(index2),EEG.pnts,EEG.trials]);
end

% command history
% -------------------
com = sprintf( '%s = pop_lms_regression( %s, [%s], %s, %s);', inputname(1), ...
    inputname(1), num2str(EOGindex), num2str(M), num2str(mu));

return;