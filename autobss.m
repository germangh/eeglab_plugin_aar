function [Y] = autobss(X, opt)
% autobss() - Performs automatic EOG artifact correction using Blind Source
% Separation (BSS) and identifying the EOG components using fractal analysis
%
% Usage:
%   >>  Y = autobss( X, opt)
%
% Inputs:
%   X              - Input data matrix, d x L
%   opt.wl         - Length of the moving-ICA analysis windows (in samples).
%                    If empty it will be set to the data length L.
%                    default: L
%   opt.ws         - Shift between correlative windows (in samples). If
%                    empty it will be set to the same value as wl.
%                    default: wl
%   opt.wl_dim     - window length (in samples) for the moving window
%                    computation of the fractal dimension
%                    default: .1*opt.wl
%   opt.ws_dim     - window shift (in samples) for the moving window
%                    computation of the fractal dimension
%                    default: opt.wl_dim
%   opt.bss_alg    - BSS algorithm to use
%                    default 'sobi'
%   opt.bss_opt    - Options to pass to the BSS algorithm
%                    default: []
%   opt.crit_alg   - criterion to use for detecting the EOG components
%                    ('fd','svf','joyce'), default: 'fd'
%   opt.crit_opt   - options to pass to the criterion function
%                    default: []
%   opt.tau        - embedding lag for computing the 'svf' criterion
%                    default: 1
%   opt.dim        - embedding dimension for computing the 'svf' criterion
%                    default: 20
%   opt.k          - index of the SVF for the 'svf' criterion
%                    default: 1
%   opt.eogindex   - indexes of the EOG channels. Necessary for computing
%                    'joyce' criterion. 
%                    default: last data channel
%
% Outputs:
%   Y   - Output data matrix (artifact corrected)
%
% Notes:
%   1) BSS will be performed on (possibly overlapping) windows of wl
%      samples.
%   2) The EOG components are detected as described in [3]. Please refer to
%      [3] when using this code in any of your publications.
%   3) Reconciling ovelapping analysis windows is done as in [1]. 
%   4) Available BSS algorithms should have an associated interface
%      function named [cmd]_ifc where [cmd] is the command name used to run
%      the BSS. A sample interface to sobi (named sobi_ifc) is included in
%      folder bss_alg.
%   5) The SVF criterion is described in [4] (not implemented yet!).
%   6) The Joyce criterion is described in [5] (not implemented yet!).
%
% References:
% [1] Wallstrom et al., G.L. International Journal of Psychology 53 (2004)
%     105-119
% [2] Katz, M., Comput. Biol. Med. 18 (1988), 145-156
% [3] Gomez-Herrero, G. , De Clercq, W., Anwar, H., Kara, O., Egiazarian, K. (2006),
%     Proceedings of NORSIG 2006, Reykjavik, Iceland.
% [4] Faul, S., Marnane, L., Lightbody, G., Boylan, G. and Connolly, S. (2005),
%     Proceedings of ICASSP 2005, Philadelphia, USA.
% [5] Joyce, C. A., Gorodnitsky, I.F., Kutas, M. (2004), Psychophysiology,
%     41, 313-325
%
%
% Author: German Gomez-Herrero <german.gomezherrero@gmail.com>
%         Institute of Signal Processing
%         Tampere University of Technology, 2007
%
% See also:
%   POP_AUTOBSSEOG, POP_AUTOBSSEMG, CMERGE_OVERLAP, CMERGE_NOOVERLAP, FD, EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com


if nargin < 1, help autoica; return; end

if ~exist('opt','var'),
    opt = def_autobss;
else
    opt = def_autobss(opt);
end

bss_alg = lower(opt.bss_alg);
ws = opt.ws;
wl = opt.wl;

% Initial and final sample of each analysis epoch
% ----------------------------------------
[d,L] = size(X);

if isempty(wl),
    wl = L; ws=L;
elseif isempty(ws),
    ws = wl;
end


% overlap between correlative samples
% -----------------------------------------
ovlength = wl-ws;

init = 1:ws:L;
final = init+wl-1;
ne = length(init);

% initialize the output (corrected)data
% ----------------------------------------
Y = zeros(d,L);

% reference EOG channels data
% ----------------------------------------
if isfield(opt.crit_opt,'eogref'),
    eogref = opt.crit_opt.eogref;
else
    eogref = [];
end

if opt.verbose,fprintf('\nRunning BSS filter in sliding windows');end
for i = 1:ne
    if final(i)>L;
        if (final(i)-L)>.5*wl,
            ne=i-1;
            %warning('Last analysis window is too short and therefore will not be processed');
            break; 
        else
            final(i)=L;
            if i < ne,
                final(i+1) = Inf;
            end
        end
    end
    Xi = X(:,init(i):final(i));
    % remove data mean
    Ximean = mean(Xi,2);
    Xi = Xi-repmat(Ximean,1,size(Xi,2));
    if ~isempty(eogref),
        eogrefi = eogref(:,init(i):final(i));
        opt.crit_opt.eogref = eogrefi;
    end

    if isfield(opt,'bss_opt') && ~isempty(opt.bss_opt),
        [Wi] = eval([bss_alg '_ifc(Xi,opt.bss_opt)']);
        Wi = real(Wi);
    else
        [Wi] = eval([bss_alg '_ifc(Xi)']);
        Wi = real(Wi);
    end
    if ~isempty(Wi),
        Ai = pinv(Wi);
        Yi = Wi*Xi;
    else
        Ai = [];
        Yi = [];
    end

    % detect the EOG-related components
    switch lower(opt.crit_alg),        
        case 'eog_fd',   % fractal dimension
            [index] = eog_fd(Yi,opt.crit_opt);
            
        case 'eog_corr', % EOG correlation
            [index] = eog_corr(Yi,opt.crit_opt);
            
        case 'emg_psd',  % PSD ratio EEGband/EMGband
            [index] = emg_psd(Yi,opt.crit_opt);

        case 'eog_svf', % singular value fraction
            error('(autobss) criterion eog_svf not working yet. Try fd instead.')
            

        case 'eog_joyce', % joyce criterion
            error('(autobss) criterion eog_joyce not working yet. Try eog_fd instead.')


    end    

    % correct the EOG (avoiding unnecessary operations and rounding errors)
    if wl>ws,
        % if there is overlap between correlative analysis windows
        if isempty(index),       % leave the data unchanged
            thisY = X(:,init(i):final(i));        
        elseif length(index)<d   % remove just some components
            EOG = real(Ai(:,index)*Yi(index,:));
            thisY = (Xi-EOG)+repmat(Ximean,1,size(Xi,2));
        else                     % remove all data
            thisY = zeros(d,wl);
        end        
        if i < 2,
            Y(:,1:final(1))=thisY;
        else
            Y(:,1:final(i)) = cmerge_overlap(Y(:,1:final(i-1)),thisY,ovlength);
        end
    else
        % if there is no overlap
        if isempty(index),       % leave the data unchanged
            Y(:,init(i):final(i)) = X(:,init(i):final(i));        
        elseif length(index)<d   % remove just some components
            EOG = real(Ai(:,index)*Yi(index,:));
            Y(:,init(i):final(i)) = (Xi-EOG)+repmat(Ximean,1,size(Xi,2));
        end                
    end

    if opt.verbose && ~mod(i,floor(ne/10)), fprintf('.'); end

end
if opt.verbose,fprintf('\nDone.\n');end

% leave the remaining data unchanged 
if final(ne)<L,
    Y(:,final(ne)+1:end) = X(:,final(ne)+1:end);
end




function [opt] = def_autobss(opt)
% def_autobss() - Sets default analysis parameters for autobss()
%
% Usage:
%   [opt] = def_autobss(opt_in)
%
% Inputs:
%   opt_in          - Input structure containing analysis options
%
% Outputs:
%   opt             - Output structure where all the analysis options not specified
%                     in opt_in have been set to their default values. The possible
%                     analysis options are:
%   opt.wl          - Length (in samples) of the analysis windows. If left
%                     empty it will be set to the data length.
%                     (default: [])
%   opt.ws          - Shift (in samples) between correlative analysis
%                     windows. If empty, the same value as wl will be
%                     assumed.
%                     (default: [])
%   opt.wl_dim      - Length (in samples) of the analysis windows used for
%                     estimating the mean fractal dimension.
%                     (default: [])
%   opt.ws_dim      - Shift (in samples) between correlative windows used for
%                     estimating the mean fractal dimension.
%                     (default: [])
%   opt.bss_alg     - BSS algorithm to use (default: 'sobi')
%   opt.bss_opt     - BSS options (default: [])
%   opt.crit_alg    - criterion to use for identifying the components to be
%                     removed. The following are available:
%                     'eog_fd','eog_corr' -> good for removing EOG
%                     'emg_psd'           -> good for removing EMG
%                     default: 'eog_fd'
%   opt.crit_opt    - Criterion options (default: [])
%   opt.tau         - embedding lag for computing the 'svf' criterion,
%                     default: 1
%   opt.dim         - embedding dimension for computing the 'svf' criterion
%                     default: 20
%   opt.k           - index of the SVF for the 'svf' criterion
%                     default: 1
%   opt.eogindex    - indexes of the EOG channels. Necessary for computing
%                     some of the criteria. 
%                     default: []
%
% Author: German Gomez-Herrero <german.gomezherrero@gmail.com>
%         Institute of Signal Processing
%         Tampere University of Technology, 2006
%
% See also:
%   POP_AUTOBSS, CMERGE_OVERLAP, CMERGE_NOOVERLAP, FD, EEGLAB
%

% Copyright (C) <2006>  <German Gomez-Herrero>
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


if ~exist('opt','var') || ~isfield(opt, 'ws'),
    opt.ws = [];
end
if ~isfield(opt, 'wl'),
    opt.wl = [];
end
if ~isfield(opt, 'wl_dim'),
    opt.wl_dim = [];
end
if ~isfield(opt, 'ws_dim'),
    opt.ws_dim = [];
end
if ~isfield(opt, 'crit_alg') || isempty(opt.crit_alg),
    opt.crit_alg = 'eog_fd';
end
if ~isfield(opt, 'crit_opt') || isempty(opt.crit_opt),
    opt.crit_opt = [];
end
if ~isfield(opt, 'tau') || isempty(opt.tau),
    opt.tau = 1;
end
if ~isfield(opt, 'dim') || isempty(opt.dim),
    opt.dim = 20;
end
if ~isfield(opt,'k') || isempty(opt.k),
    opt.k = 1;
end
if ~isfield(opt, 'verbose') || isempty(opt.verbose),
    opt.verbose = 1;
end
if ~isfield(opt,'bss_alg') || isempty(opt.bss_alg),
    opt.bss_alg = 'sobi';
end
if ~isfield(opt,'bss_opt') || isempty(opt.bss_opt),
    opt.bss_opt = [];
end
