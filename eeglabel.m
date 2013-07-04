function [EEG,com] = eeglabel(EEG,regions)
% eeglabel() - labels portions of continuous data in an EEGLAB dataset
%
% Usage:
%   >> EEGOUT = eeglabel(EEGIN, regions)
%
% Inputs:
%   INEEG      - input dataset
%   regions    - array of regions to suppress. number x [beg end]  of
%                regions. 'beg' and 'end' are expressed in term of points
%                in the input dataset. Size of the array is
%                number x 2 of regions.
%
% Outputs:
%   INEEG      - output dataset with updated data, events latencies and
%                additional events.
%
% Author: German Gomez-Herrero <german.gomezherrero@tut.fi>
%         Institute of Signal Processing
%         Tampere University of Technology, 2008
%
% See also:
%   POP_EEGLABEL, EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

if nargin < 2,
    help eeglabel;
end
if isempty(regions),
    return;
end

% open a window to get the label value
% --------------------------------------
uigeom = {[1.5 1];[1.5 1]};
uilist = {{'style' 'text' 'string' 'Label for this EEG epoch(s):'} ...
    {'style' 'edit' 'string' ''} ...
    {'style' 'text' 'string' 'Short description:'} ...
    {'style' 'edit' 'string' ''} ...
    };
guititle = 'Choose a label - eeglabel()';
result = inputgui( uigeom, uilist, 'pophelp(''eeglabel'')', guititle, [], 'normal');

label = eval(['''' result{1} '''']);
descr = eval(['''' result{2} '''']);


% handle regions from eegplot and insert labels
% -------------------------------------
if size(regions,2) > 2,
    regions = regions(:,3:4);
end
for i = 1:size(regions,1),
    startevent.type = [label '.start'];
    endevent.type = [label '.end'];
    startevent.latency = regions(i,1)-.5;
    endevent.latency = regions(i,2)-.5;
    startevent.duration = regions(i,2)-regions(i,1);
    endevent.duration = 5;
    startevent.description = descr;
    endevent.description = descr;    
    EEG.event = [EEG.event startevent endevent];
end
    
EEG = eeg_checkset(EEG,'eventconsistency');

com = sprintf('%s = eeglabel( %s, %s);', inputname(1), inputname(1), vararg2str({ regions }));
return;
