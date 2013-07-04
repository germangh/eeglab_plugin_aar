function [index,I] = eog_fd(X,opt)
% eog_fd() - Selects EOG components according to their fractal dimensions
%
% Usage:
%   >> [index] = eog_fd(X,opt)
%
% Inputs:
%   X           - data matrix (dxN, data channels are rowwise)
%   opt.wl      - window length for computing the mean fractal dimension
%                 def: floor(.1*N)
%   opt.ws      - window shift for computing the mean fractal dimension
%                 def: wl
%   opt.method  - method for computing the fractal dimension
%                 ('sevcik','katz','sevcik_mean','sevcik_var','katz_mean',
%                 'katz_var')
%                 def: 'sevcik_mean';
%   opt.range   - range of components that can be removed. At least
%                 opt.range(1) components will be removed in each analysis
%                 window and at most opt.range(2) components.
%                 def: [min(2,floor(d/9)) floor(d/3)]
%
% Outputs:
%   index   - indexes of the rows of X corresponding to EOG components
%
% References:
% [1] G. Gomez-Herrero, W. De Clercq, H. Anwar, O. Kara, K. Egiazarian,
% S. Van Huffel, W. Van Paesschen. Automatic removal of ocular artifacts in
% the EEG without an EOG reference channel, Proceedings of NORSIG 2006.
%
%
% See also:
%   FD, EOG_CORR, POP_AUTOBSSEOG, AUTOBSS, EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com


if nargin < 1, help eog_fd; return; end
[d,N] = size(X);
if ~exist('opt','var'),
    opt = def_eog_fd;
else
    opt = def_eog_fd(opt);
end

% default range of components that might be removed
% ---------------------------------------------
if isempty(opt.range),
    RANGE = [min(2,floor(d/9)) floor(d/3)]; 
else
    RANGE = opt.range;
end


method = opt.method;
if isempty(opt.wl),
    wl = floor(.1*N);
end
if isempty(opt.ws),
    ws = wl;
end

% compute fractal dimension of components
% ---------------------------------------------
fractal_dim = zeros(1,d);
fdopt.method = method;
fdopt.wl = wl;
fdopt.ws = ws;
for j = 1:d    
    fractal_dim(j) = fd(X(j,:),fdopt);
end

% sort components by increasing fractal dimension
% ---------------------------------------------
[fractal_dim,I] = sort(fractal_dim);
fdistance = fractal_dim(2:end)-fractal_dim(1:end-1);

% separe EOG components from non-EOG components
% ---------------------------------------------
index = I(1);
for j = 1:(length(fdistance)-1)
    if mean(cumsum(fdistance(j+1:end)))<mean(cumsum(fdistance(1:j))),
        index = I(1:j);
        break;
    end
end

% take a number of components within the specified range
% ---------------------------------------------
if length(index) < RANGE(1), index = I(1:RANGE(1)); end
if length(index) > RANGE(2), index = I(1:RANGE(2)); end

return;



% subfunction to define the default parameters
% --------------------------------------------
function [opt] = def_eog_fd(opt)

if nargin < 1 || isempty(opt) || ~isfield(opt,'wl'),
    opt.wl = [];
end
if ~isfield(opt,'ws'),
    opt.ws = [];
end
if ~isfield(opt,'method') || isempty(opt.method),
    opt.method = 'sevcik_mean';
end
if ~isfield(opt,'range'),
    opt.range = [];
end

