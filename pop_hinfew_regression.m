% pop_hinftv_regression() - Automatic EOG correction using H infinity norm
% time-varying algorithm described in [1]
%
% Usage:
%   >> [OUTEEG] = pop_hinftv_regression(INEEG, EOGchans, M, eta, rho, eps, lambda, evchans)
%
% Inputs:
%   INEEG    - input EEG dataset
%   EOGchans - indexes of the EOG reference channels
%   M        - Order of the adaptive filter
%   eta      - Factor reflecting a priori knowledge of how close the
%              estimated filter weights at t=0 are to their optimal
%              value at that time instant (eta << 1)
%   rho      - Factor reflecting a priori knowledge of how rapidly
%              the filter coefficients vary with time
%   eps      - Positive constant described in [1] (~ 1)
%   lambda   - Forgetting factor
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
%
% Author: German Gomez-Herrero
%         <german.gomezherrero@ieee.org>
%         Institute of Signal Processing
%         Tampere University of Technology, 2007
%
% References:
% [1] S. Puthusserypady and T. Ratnarajah, IEEE Signal Processing Letters
% 12, 816-819
%
% See also:
%   HINFEW_REGRESSION, POP_HINFTV_REGRESSION, HINFTV_REGRESSION, EEGLAB

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

function [EEG,com] = pop_hinfew_regression(EEG, EOGindex, M, eta, rho, eps, lambda, evchans)

com = '';

if nargin < 1,
    help pop_hinfew_regression;
    return;
end
if nargin < 2, EOGindex = []; end
if nargin < 3, M = 3; end
if nargin < 4, eta = 5e-3; end
if nargin < 5, rho = 1e-5; end
if nargin < 6, eps = 1.5; end
if nargin < 7, lambda = 0.99; end
if nargin < 8, evchans = []; end

if isempty(EEG.data)
    disp('(pop_hinfew_regression) error: cannot clean an empty dataset'); 
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
if nargin < 8,    
    uigeom = {[1.5 1] [1.5 1] [1.5 1] [1.5 1] [1.5 1] [1.5 1] [1.5 1]};
    uilist = { { 'style' 'text'      'string'    'EOG channel indexes:'} ...
        {'style'  'edit' 'string'    num2str(EOGindex)} ...
        {'style'  'text'      'string'    'Filter order (M):'} ...
        {'style'  'edit'      'string'    num2str(M)} ...
        {'style'  'text'      'string'    'Distance at t=0 to optimal solution (eta):'} ...
        {'style'  'edit'      'string'    num2str(eta)} ...        
        {'style'  'text'      'string'     'Speed of variation of filter coefficients (rho):'} ...
        {'style'  'edit'      'string'    num2str(rho)} ...
        {'style'  'text'      'string'     'Positive constant epsilon:'} ...
        {'style'  'edit'      'string'    num2str(eps)} ...
        {'style'  'text'      'string'     'Forgetting factor (lambda):'} ...
        {'style'  'edit'      'string'    num2str(lambda)} ...
        {'style'  'text'      'string'     'Store filter weights for channels:'} ...
        {'style'  'edit'      'string'    num2str(evchans)} ...  
        };
    guititle = 'Correct EOG using HinfEW regression -- pop_hinfew_regression()';
    result = inputgui( uigeom, uilist, 'pophelp(''pop_hinfew_regression'')', guititle, [], 'normal');

    if isempty(result), return; end;

    % reading params
    % -------------------
    EOGindex = eval(['[' result{1} ']']);
    M = eval(['[' result{2} ']']);
    eta = eval(['[' result{3} ']']);    
    rho = eval(['[' result{4} ']']);
    eps = eval(['[' result{5} ']']);
    lambda = eval(['[' result{6} ']']);
    evchans = eval(['[' result{7} ']']);
end

% build the options structure
% -----------------------------
opt.refdata = reshape(EEG.data(EOGindex,:,:),length(EOGindex),EEG.pnts*EEG.trials);
opt.M = M;
opt.eta = eta;
opt.rho = rho;
opt.eps = eps;
opt.lambda = lambda;
EEGindex = find(~ismember(1:EEG.nbchan,EOGindex));

% run the EOG correction
% ----------------------
tmp_in = reshape(EEG.data,EEG.nbchan,EEG.pnts*EEG.trials);
index1 = intersect(EEGindex,evchans);
index2 = intersect(EEGindex,find(~ismember(1:EEG.nbchan,evchans)));
if ~isempty(index1),
    [tmp,H,Hh] = hinfew_regression(tmp_in(index1,:), opt);
    EEG.data(index1,:,:) = reshape(tmp,[length(index1),EEG.pnts,EEG.trials]);
    EEG.Hh = Hh;
    EEG.evchans = evchans;
end
if ~isempty(index2),
    tmp = hinfew_regression(tmp_in(index2,:), opt);
    EEG.data(index2,:,:) = reshape(tmp,[length(index2),EEG.pnts,EEG.trials]);
end


% command history
% -------------------
com = sprintf( '%s = pop_hinfew_regression( %s, [%s], %s, %s, %s, %s, %s, [%s]);', ...
    inputname(1), inputname(1), num2str(EOGindex), num2str(M), num2str(eta), ...
    num2str(rho), num2str(eps), num2str(lambda),num2str(evchans));

return;