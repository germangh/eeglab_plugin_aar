function [index,I] = eog_corr(X,opt)
% eog_corrcoef() - Selects EOG components according to their correlation
%                  with the EOG channels
%
% Usage:
%   >> [index] = eog_corr(X,opt)
%
% Inputs:
%   X               - data matrix (dxN, data channels are rowwise)
%   opt.eogref      - EOG data matrix (kxN, data channels are rowwise)
%   opt.corrcoef    - Absolute correlation between a component and any of
%                     the EOG channels to label the component as artifactual
%                     (default: 0.75)
%   opt.range       - range of components that can be removed. At least
%                     opt.range(1) components will be removed in each analysis
%                     window and at most opt.range(2) components.
%                     def: [min(2,floor(d/9)) floor(d/3)]
%                    
% Outputs:
%   index   - indexes of the rows of X corresponding to EOG components
%
% See also:
%   EOG_FD, POP_AUTOBSSEOG, AUTOBSS, EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

if nargin < 1, help eog_corr; return; end
[d] = size(X,1);
if ~exist('opt','var'),
    opt = def_eog_corr;
else
    opt = def_eog_corr(opt);
end

% default range of components that might be removed
% ---------------------------------------------
if isempty(opt.range),
    RANGE = [min(2,floor(d/9)) floor(d/3)]; 
else
    RANGE = opt.range;
end

% EOG reference channels
% ---------------------------------------------
eogref = opt.eogref;

% compute maximum EOG correlation of each component
% ---------------------------------------------
R = zeros(d,1);
for j = 1:d    
    for k = 1:size(eogref,1)       
        R(j) = max(R(j),corr2(squeeze(X(j,:)),squeeze(eogref(k,:))).^2);        
    end
end

% sort by increasing correlation with EOG channel(s)
% ---------------------------------------------
[R,I] = sort(R,'descend');
index = find(R>opt.corrcoef);
if isempty(index), index = I(1); end

% take a number of components within the specified range
% ---------------------------------------------
if length(index) < RANGE(1), index = I(1:RANGE(1)); end
if length(index) > RANGE(2), index = I(1:RANGE(2)); end

return;

% subfunction to define the default parameters
% --------------------------------------------
function [opt] = def_eog_corr(opt)

if ~exist('opt','var') || ~isfield(opt, 'eogref') || isempty(opt.eogref),
    error('(eog_corr) EOG reference data must be provided')
end
if ~isfield(opt, 'verbose') || isempty(opt.verbose),
    opt.verbose = 1;
end
if ~isfield(opt,'corrcoef') || isempty(opt.corrcoef),
    opt.corrcoef = 0.75;
end
if ~isfield(opt,'range'),
    opt.range = [];
end

