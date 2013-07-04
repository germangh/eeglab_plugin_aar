function [index,I] = emg_psd(X,opt)
% emg_psd() - Selects EMG components according to their Power Density
% Spectrum
%
% Usage:
%   >> [index,I] = emg_psd(X,opt)
%
% Inputs:
%   X               - data matrix (dxN, data channels are rowwise)
%   opt.femg       - frequency aproximately separating the EEG and EMG bands
%                    def: 15 Hz
%   opt.fs         - sampling frequency
%                    def: 250 Hz
%   opt.ratio      - ratio of average power (per unit of frequency) in EEG
%                    band to average power in the EMG band below which a
%                    component will be considered to be EMG-related.
%                    Increasing/decreasing this value increases/decreases
%                    the amount of correction.
%                    def: 10
%   opt.range      - range of components that can be removed. At least
%                    opt.range(1) components will be removed in each analysis
%                    window and at most opt.range(2) components.
%                    def: [0 floor(d/2)]
%   opt.estimator  - espectral estimator object. This can be any of the
%                    spectral estimatio objects provided by the constructor
%                    spectrum.estmethod. By default a spectrum.welch
%                    estimator is used with segment length larger or equal
%                    than half the data points of a component and larger or
%                    equal than 2 seconds (the minimum of these two values
%                    will be used as segment length). If opt.estimator is
%                    left empty or the field .estimator does not exist in
%                    the opt struct then the default Welch estimator will be
%                    used.
%   opt.NFFT       - specifies the number of FFT points to use to calculate
%                    power espectral density. By default this is set to be
%                    equal to the segment length of the default spectral
%                    estimator (see above). If left empty or if the field
%                    .NFFT does not exist then the default value will be
%                    used.
%
% Outputs:
%   index   - indexes of the rows of X corresponding to EMG components
%
%
% Notes:
%   - This function requires MATLAB Signal Processing Toolbox (SPT) v. 6.2
%     or newer, i.e. the version of the SPT included with MATLAB v. 7.0
%     (R14) or newer MATLAB releases.
%   - The optimum value of the parameter opt.ratio depends on the actual
%     spectral estimator used.
%
% See also:
%   FD, POP_AUTOBSSEMG, AUTOBSS, EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com


% check version of the Signal Processing Toolbox
% -------------------------------------------------
sptver = ver('signal');
if isempty(sptver),
    error('(emg_psd) this criterion requires Signal Processing Toolbox');
elseif (str2double(sptver.Version)<6.2 && datenum(sptver.Date)<datenum('1-Jan-2007')),
    error('(emg_psd) this criterion needs the Signal Processing Toolbox v.6.2 or newer');
elseif (str2double(sptver.Version)< 6.6 && datenum(sptver.Date)<datenum('1-Jan-2007')),
    sptold = true;
else
    sptold = false;
end

if nargin < 1, help emg_psd; return; end
[d,N] = size(X);

% remove mean from data
% ------------------------------------------------
X = X-repmat(mean(X,2),1,N);

if ~exist('opt','var'),
    opt = def_emg_psd;
else
    opt = def_emg_psd(opt);
end
fs       = opt.fs;
ratio_th = opt.ratio;
femg     = opt.femg;

% default range of components that might be removed
% -------------------------------------------------
if isempty(opt.range),
    RANGE = [0 floor(d/2)];
else
    RANGE = opt.range;
end

% default spectral estimator
% -------------------------------------------------
if isempty(opt.estimator),
    sl = min(2*opt.fs,floor(N/2));
    if sptold,
        h = spectrum.welch({'Hamming'},sl);        
    else
        h = spectrum.welch({'Hamming'},sl);
    end
else
    h=opt.estimator;
end

% default FFT length
% -------------------------------------------------
if isempty(opt.NFFT),
    NFFT = 2^ceil(log2(min(2*opt.fs,floor(N/2))));
else
    NFFT = opt.NFFT;
end

if sptold && ~isempty(NFFT),
    set(h,'FFTLength','UserDefined');
end

% Compute the average power below and over femg
% -------------------------------------------------
p1 = zeros(1,d);
p2 = zeros(1,d);
for i = 1:d
    hpsd = psd(h,X(i,:),'NFFT',NFFT);
    p1(i) = avgpower(hpsd,[0 femg/(fs/2)*pi])/(femg/(fs/2)*pi);
    p2(i) = avgpower(hpsd,[femg/(fs/2)*pi pi])/(pi-femg/(fs/2)*pi);
end

% detect components that are likely to be EMG-related
% ---------------------------------------------------
ratio = p1./p2;
[ratio,I] = sort(ratio);
index=I(ratio<ratio_th);


% take a number of components within the specified range
% ---------------------------------------------
if length(index) < RANGE(1), index = I(1:RANGE(1)); end
if length(index) > RANGE(2), index = I(1:RANGE(2)); end

return;


% subfunction to define the default parameters
% ---------------------------------------------------
function [opt] = def_emg_psd(opt)
if nargin < 1 || ~isfield(opt,'fs'),
    opt.fs=250;
    warning('(emg_psd) Default sampling frequency [%d] Hz will be used',opt.fs);
end
if ~isfield(opt, 'bssout'),
    opt.bssout = [];
end
if ~isfield(opt, 'ratio') || isempty(opt.ratio),
    opt.ratio = 10;
end
if ~isfield(opt, 'femg') || isempty(opt.femg),
    opt.femg = 15;
end
if ~isfield(opt,'range'),
    opt.range = [];
end
if ~isfield(opt,'estimator'),
    opt.estimator = [];
end
if ~isfield(opt,'NFFT'),
    opt.NFFT = [];
end

