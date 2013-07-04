% pop_autobssemg() - Automatic EMG correction using Blind Source Separation
%
% Usage:
%   >> OUTEEG = pop_autobssemg(INEEG, wl, ws, bss_alg, bss_opt, crit_alg, crit_opt)
%
% Inputs:
%   INEEG    - input EEG dataset
%   wl       - analysis window length (in seconds)
%   ws       - shift between correlative analysis windows (in seconds)
%   bss_alg  - name of the BSS algorithm to use
%   bss_opt  - options to pass to the BSS algorithm as
%              {opt_name1, opt_value1,opt_name2,opt_value2,...}
%   crit_alg - name of the criterion for selecting the components to remove
%   crit_opt  - options to pass to the criterion function as
%              {opt_name1, opt_value1,opt_name2,opt_value2,...}
%
% Outputs:
%   OUTEEG  - output dataset
%
% Notes:
%
% (*) The automatic method for correcting EMG artifacts using CCA was
%     originally proposed by Wim et al. in:
%
%     De Clercq, W. et al., A new muscle artifact removal technique to improve
%     interpretation of the ictal scalp electroencephalogram, Proceedings of
%     EMBC 2005, Shanghai, China, pp. 1136-1139.
%
% (*) When using pop_autobssemg from the command line a pop_up window will
%     be displayed ONLY if the user did not provide the value of the window
%     length, that is, only if there is a single input parameter. 
%
% See also:
%   AUTOBSS, EMG_PSD, EEGLAB


% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

function [EEG, com] = pop_autobssemg(EEG, wl, ws, bss_alg, bss_opt, crit_alg, crit_opt)

com = '';
if nargin < 2,
    % at least the first 2 params are required or the pop-up window will appear
    showpopup = true;
else
    showpopup = false; % do not show pop-up window unless completely necessary
end

if nargin < 1,
    help pop_autobssemg;
    return;
end

if isempty(EEG.data)
    disp('(pop_autobssemg) error: cannot clean an empty dataset'); return;
end;

% default window length/shift (just a rule of thumb)
% ---------------------------
N = EEG.pnts*size(EEG.data,3);
def_wl_samples = EEG.srate*0.02*EEG.nbchan^2;
if size(EEG.data,3) > 1,
    def_wl_samples = def_wl_samples-mod(def_wl_samples,EEG.pnts);
end
def_wl_samples = min(def_wl_samples,N);
def_wl = def_wl_samples/EEG.srate;
def_ws = def_wl;

% use default wl and ws?
% ---------------------------
if nargin < 2 || (isempty(wl) && isempty(ws)),
    wl = def_wl;
    ws = def_ws;
elseif nargin < 3,
    ws = wl;
end

% find available algorithms
% -----------------------------
allalgs = {'bsscca','iwasobi','efica','multicombi','fcombi','sobi','fpica', ...
    'runica','jader','fastica','pca'};
selectalg = {};

for index = 1:length(allalgs)
    if exist([allalgs{index} '_ifc'],'file') && exist([allalgs{index}],'file'),
        selectalg = {selectalg{:} allalgs{index}};
    end
end
if isempty(selectalg),
    error('(pop_autobssemg) I could not find an interface function for any BSS algorithm');
end

% use default BSS algorithm with default options?
% ----------------------------------------------------
if nargin < 4 || isempty(bss_alg),
    bss_alg = selectalg{1};
elseif isempty(bss_alg) && nargin > 4 && ~isempty(bss_opt),
    showpopup = true; % wrong input -> show pop-up window!
end
if nargin < 5 || isempty(bss_opt),  
   switch lower(bss_alg),
       case 'bsscca',
           bss_opt = {'eigratio',1e6};           
       otherwise,
           bss_opt = {};           
   end
end

% find available criteria to remove components
% ---------------------------------------------
sptver = ver('signal');

if isempty(sptver) || (str2double(sptver.Version)<6.2 && datenum(sptver.Date)<datenum('1-Jan-2007')),
    error('(pop_autobssemg) all available criteria require MATLAB''s Signal Processing Toolbox v.6.2 or newer');
else
    allcrits = {'emg_psd'};
end
selectcrit = {};
for index = 1:length(allcrits)
    if exist([allcrits{index}],'file') && exist([allcrits{index}],'file'),
        selectcrit = {selectcrit{:} allcrits{index}};
    end
end
if isempty(selectcrit),
    error('(pop_autobssemg) I could not find any criterion function');
end

% user wants to use the default criterion with default options
% ----------------------------------------------------
if nargin < 6 || isempty(crit_alg),
    crit_alg = selectcrit{1};
elseif isempty(crit_alg) && nargin > 6 && ~isempty(crit_alg),
    showpopup = true; % wrong input -> show pop-up window!
end

if nargin < 7 || isempty(crit_opt),
    % user wants def. criterion options    
    switch lower(crit_alg),
        case 'emg_psd',
            sl = min(floor(EEG.srate/2),floor(wl/5)*EEG.srate); % 5 windows will be used by default
            estimator = spectrum.welch({'Hamming'},sl);
            crit_opt = {'ratio',10,'fs',num2str(EEG.srate),'femg',15,'estimator',...
                estimator,'range',[0 floor(EEG.nbchan/2)]};
        otherwise,
            crit_opt = {};
    end
end

% default spectral estimator (assumes default criterion is emg_psd)
% ---------------------------
if ~strcmpi(selectcrit{1},'emg_spd'),
    sl = min(floor(EEG.srate/2),floor(wl/5)*EEG.srate); % 5 windows will be used by default
    def_estimator = ['spectrum.welch({''Hamming''},' num2str(sl) ')'];           
else
    def_estimator = '';
end

% build the string of default criterion options
% ---------------------------
switch lower(selectcrit{1}),
    case 'emg_psd',
        def_criterion_opts = ['''ratio'',10,''fs'',' num2str(EEG.srate) ...
    ',''femg'',15,''estimator'',' def_estimator ',''range'',[' ... 
        num2str([0 floor(EEG.nbchan/2)]) ']'];
    otherwise,
        def_criterion_opts = '';        
end


% display input dialog
% ---------------------------
if showpopup,
    uigeom = {[1.5 1] [1.5 1] [1.5 1] 1 1 1 [1.5 1] 1 1};
    uilist = { { 'style' 'text'      'string'    'BSS algorithm:'} ...
        { 'style' 'popupmenu' 'string'    selectalg} ...
        {'style'  'text'      'string'    'Analysis window length (seconds):'} ...
        {'style'  'edit'      'string'    num2str(def_wl)} ...
        {'style'  'text'      'string'    'Shift between correlative windows (seconds):'} ...
        {'style'  'edit'      'string'    num2str(def_ws)} ...
        {'style'  'text'      'string'    'Options to pass to the BSS algorithm ([option_name],[value],...):'} ...
        {'style'  'edit'      'string'    '''eigratio'',1e6'} ...
        {'style'  'text'      'string'    ''} ...
        { 'style' 'text'      'string'    'Criterion to remove components:'} ...
        { 'style' 'popupmenu' 'string'    selectcrit} ...
        {'style'  'text'      'string'    'Options to pass to the criterion ([option_name],[value],...):'} ...
        {'style'  'edit'      'string'    def_criterion_opts} ...
        };
    guititle = 'Correct EMG using BSS -- pop_autobssemg()';
    result = inputgui( uigeom, uilist, 'pophelp(''pop_autobssemg'')', guititle, [], 'normal');

    if isempty(result), return; end;

    % reading params
    % -------------------
    bss_alg = selectalg{result{1}};
    wl = eval(['[' result{2} ']']);
    ws = eval(['[' result{3} ']']);    

    % read BSS parameters
    % -------------------
    bss_opt = eval(['{' result{4} '}']);

    % criterion
    % -----------
    crit_alg = selectcrit{result{5}};

    % criterion parameters
    % ---------------------
    crit_opt = eval(['{' result{6} '}']);   
end

% correct wl,ws so that they will be an integer number of trials
% -------------------
wl = min(floor(wl*EEG.srate),EEG.pnts*EEG.trials);
ws = min(floor(ws*EEG.srate),EEG.pnts*EEG.trials);
if size(EEG.data,3) > 1,
    wl = wl-mod(wl,EEG.pnts);
    ws = ws-mod(ws,EEG.pnts);
end

opt = struct;

% build BSS parameters structure
% -------------------------------
if exist('bss_opt','var'),
    opt.bss_opt = struct;
    for i = 1:2:length(bss_opt),
        opt.bss_opt.(bss_opt{i}) = bss_opt{i+1};
    end
end

% build criterion parameters structure
% -------------------------------
if exist('crit_opt','var'),
    opt.crit_opt = struct;
    for i = 1:2:length(crit_opt),
        opt.crit_opt.(crit_opt{i}) = crit_opt{i+1};
    end
end
opt.crit_opt.eogref = [];

% rest of parameters
% -------------------
opt.wl = wl;
opt.ws = ws;
opt.bss_alg = bss_alg;
opt.crit_alg = crit_alg;
opt.crit_opt.fs = EEG.srate;

% run the EMG correction
% -------------------
EEG.data = autobss(reshape(EEG.data,EEG.nbchan,EEG.pnts*EEG.trials),opt);
EEG.data = reshape(EEG.data,[EEG.nbchan,EEG.pnts,EEG.trials]);

% command history
% -------------------
bss_opt_str = [];
for i = 1:2:length(bss_opt),
    if isnumeric(bss_opt{i+1}),
        bss_opt_str = [bss_opt_str ',''' bss_opt{i} ''', [' num2str(bss_opt{i+1}) ']'];
    elseif ischar(bss_opt{i+1}),
        bss_opt_str = [bss_opt_str ',''' bss_opt{i} ''',''' bss_opt{i+1} ''''];
    else
        bss_opt_str = [bss_opt_str ',''' bss_opt{i} ''',''' class(bss_opt{i+1}) ''''];        
    end
end
if ~isempty(bss_opt_str),
    bss_opt_str(1)=[];
end

crit_opt_str = [];
for i = 1:2:length(crit_opt),
    if isnumeric(crit_opt{i+1}),
        crit_opt_str = [crit_opt_str ',''' crit_opt{i} ''', [' num2str(crit_opt{i+1}) ']'];
    elseif ischar(crit_opt{i+1}),
        crit_opt_str = [crit_opt_str ',''' crit_opt{i} ''',''' crit_opt{i+1} ''''];
    else
        crit_opt_str = [crit_opt_str ',''' crit_opt{i} ''',' class(crit_opt{i+1})];
    end
end
if ~isempty(crit_opt_str),
    crit_opt_str(1)=[];
end

com = sprintf( '%s = pop_autobssemg( %s, [%s], [%s], ''%s'', {%s}, ''%s'', {%s});', inputname(1), ...
    inputname(1), num2str(wl/EEG.srate), num2str(ws/EEG.srate), bss_alg, ...
    bss_opt_str, crit_alg, crit_opt_str);

return;
