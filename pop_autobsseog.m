% pop_autobsseog() - Automatic EOG correction using Blind Source Separation
%
% Usage:
%   >> OUTEEG = pop_autobsseog(INEEG, wl, ws, bss_alg, bss_opt, crit_alg, crit_opt)
%
% Inputs:
%   INEEG    - input EEG dataset
%   wl       - analysis window length
%   ws       - shift between correlative analysis windows
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
% (*) The automatic method for correcting EOG artifacts using BSS and
%     and criterion 'eog_fd' (based in fractal analysis) was first proposed
%     in:
%
%     Gómez-Herrero, G. et al., Automatic removal of ocular artifacts in
%     the EEG without a reference EOG channel, in Proc. NORSIG 2006,
%     Reykjavik, Iceland, pp. 130-133, 2006.
%
% (*) When using pop_autobssemg from the command line a pop_up window will
%     be displayed ONLY if the user did not provide the value of the window
%     length, that is, only if there is a single input parameter. 
%
%
% See also:
%   AUTOBSS, EOG_FD, EOG_CORR, EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

function [EEG, com] = pop_autobsseog(EEG, wl, ws, bss_alg, bss_opt, crit_alg, crit_opt)

com = '';
if nargin < 2,
    % at least the first 2 params are required or the pop-up window will appear
    showpopup = true; 
else 
    showpopup = false; % do not show pop-up window unless necessary
end

if nargin < 1,
    help pop_autobsseog;
    return;
end
if isempty(EEG.data)
    disp('(pop_autobsseog) error: cannot clean an empty dataset');
    return;
end;

% default window length/shift (just a rule of thumb)
% ---------------------------
N = EEG.pnts*size(EEG.data,3);
def_wl_samples = EEG.srate*0.5*EEG.nbchan^2;
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
allalgs = {'sobi','iwasobi','efica','multicombi','fcombi','fpica', ...
    'runica','jader','fastica','bsscca','pca'};
selectalg = {};
for index = 1:length(allalgs)
    if exist([allalgs{index} '_ifc'],'file') && exist([allalgs{index}],'file'),
        selectalg = {selectalg{:} allalgs{index}};
    end
end
if isempty(selectalg),
    error('(pop_autobsseog) I could not find an interface function for any BSS algorithm');
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
       case 'sobi',
           bss_opt = {'eigratio',1e6};           
       otherwise,
           bss_opt = {};           
   end
end

% find available criteria to remove components
% ---------------------------------------------
allcrits = {'eog_fd','eog_corr','eog_joyce','eog_svf'};
selectcrit = {};
for index = 1:length(allcrits)
    if exist([allcrits{index}],'file'),
        selectcrit = {selectcrit{:} allcrits{index}};
    end
end
if isempty(selectcrit),
    error('(pop_autobsseog) I could not find any criterion function');
end

% default minimum/maximum number of components to be removed
% -----------------------------------------------------------
def_mincmp = min(2,floor(EEG.nbchan/9));
def_maxcmp = max(def_mincmp,floor(EEG.nbchan/3));

% use default criterion with default options?
% ----------------------------------------------------
if nargin < 6 || isempty(crit_alg),
    crit_alg = selectcrit{1};
elseif isempty(crit_alg) && nargin > 6 && ~isempty(crit_alg),
    showpopup = true; % wrong input -> show pop-up window!
end

if nargin < 7 || isempty(crit_opt),
    % user wants def. criterion options    
    switch lower(crit_alg),
        case 'eog_fd',            
            crit_opt = {'range',[def_mincmp def_maxcmp]};
        otherwise,
            crit_opt = {};
    end
end


% try to guess which are the EOG channels (if any)
% -------------------------------------------------
EOGindex = [];
for i = 1:length(EEG.chanlocs),
    labels = EEG.chanlocs(i).labels;
    if ~isempty(strfind(lower(labels),'eog')),
        EOGindex = [EOGindex i];
    end
end



% display input dialog
% ---------------------------
if showpopup,
    uigeom = {[1.5 1] [1.5 1] [1.5 1] 1 1 1 [1.5 1] [1.5 1] 1 1};
    uilist = { { 'style' 'text'      'string'    'BSS algorithm:'} ...
        { 'style' 'popupmenu' 'string'    selectalg} ...
        {'style'  'text'      'string'    'Analysis window length (seconds):'} ...
        {'style'  'edit'      'string'    num2str(def_wl)} ...
        {'style'  'text'      'string'    'Shift between correlative windows (seconds):'} ...
        {'style'  'edit'      'string'    num2str(def_ws)} ...
        {'style'  'text'      'string'    'Options to pass to the BSS algorithm ([option_name],[value],...):'} ...
        {'style'  'edit'      'string'    '''eigratio'',1e6'} ...
        {'style'  'text'      'string'    ''} ...
        {'style'  'text'      'string'    'Criterion to remove components:'} ...
        {'style' 'popupmenu' 'string'    selectcrit} ...
        {'style'  'text'      'string'    'EOG channels indexes:'} ...
        {'style'  'edit'      'string'    num2str(EOGindex)} ...
        {'style'  'text'      'string'    'Options to pass to the criterion ([option_name],[value],...):'} ...
        {'style'  'edit'      'string'    ['''range'',[' num2str(def_mincmp) ',' num2str(def_maxcmp) ']']} ...
        };
    guititle = 'Correct EOG using BSS -- pop_autobsseog()';
    result = inputgui( uigeom, uilist, 'pophelp(''pop_autobsseog'')', guititle, [], 'normal');

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

    % EOG channels indexes
    % ---------------------
    EOGindex = eval(['[' result{6} ']']);

    % criterion parameters
    % ---------------------
    crit_opt = eval(['{' result{7} '}']);

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
        opt.bss_opt.(bss_opt{i})= bss_opt{i+1};
    end
end

% build criterion parameters structure
% -------------------------------
if exist('crit_opt','var'),
    opt.crit_opt = struct;
    for i = 1:2:length(crit_opt),
        opt.crit_opt.(crit_opt{i}) =  crit_opt{i+1};
    end
end
opt.crit_opt.eogref = reshape(EEG.data(EOGindex,:),length(EOGindex),EEG.pnts*EEG.trials);

% rest of parameters
% -------------------
opt.wl = wl;
opt.ws = ws;
opt.bss_alg = bss_alg;
opt.crit_alg = crit_alg;

% run the EOG correction
% -------------------
EEG.data = autobss(reshape(EEG.data,EEG.nbchan,EEG.pnts*EEG.trials),opt);
EEG.data = reshape(EEG.data,[EEG.nbchan,EEG.pnts,EEG.trials]);

% command history
% -------------------
bss_opt_str = [];
for i = 1:2:length(bss_opt),
    if isnumeric(bss_opt{i+1}),
        bss_opt_str = [bss_opt_str ',''' bss_opt{i} ''', [' num2str(bss_opt{i+1}) ']'];
    else
        bss_opt_str = [bss_opt_str ',''' bss_opt{i} ''',''' bss_opt{i+1} ''''];
    end
end
if ~isempty(bss_opt_str),
    bss_opt_str(1)=[];
end

crit_opt_str = [];
for i = 1:2:length(crit_opt),
    if isnumeric(crit_opt{i+1}),
        crit_opt_str = [crit_opt_str ',''' crit_opt{i} ''',[' num2str(crit_opt{i+1}) ']'];
    else
        crit_opt_str = [crit_opt_str ',''' crit_opt{i} ''',''' crit_opt{i+1} ''''];
    end
end
if ~isempty(crit_opt_str),
    crit_opt_str(1)=[];
end

com = sprintf( '%s = pop_autobsseog( %s, [%s], [%s], ''%s'', {%s}, ''%s'', {%s});', inputname(1), ...
    inputname(1), num2str(wl/EEG.srate), num2str(ws/EEG.srate), bss_alg, ...
    bss_opt_str, crit_alg, crit_opt_str);

return;
