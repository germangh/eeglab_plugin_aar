function [com] = pop_eeglabel(EEG)
% pop_eeglabel() -  Opens an EEGPLOT window with a button that allowes the
% labeling of EEG periods. 
%
% Usage:
%   >> [com] = pop_eeglabel(EEG)
%
% Inputs:
%   EEG - EEGLAB dataset structure
%
% Outputs:
%   com - The equivalent command line command
%
% Author: German Gomez-Herrero <german.gomezherrero@tut.fi>
%         Institute of Signal Processing
%         Tampere University of Technology, 2008
%
% See also:
%   EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com


com = '';
if nargin < 1, 
    help pop_eeglabel;
    return;
end

command = ...
    [ '[EEGTMP LASTCOM] = eeglabel(EEG,eegplot2event(TMPREJ,-1));',...
    'if ~isempty(LASTCOM),' ...
    '  [ALLEEG EEG CURRENTSET tmpcom] = pop_newset(ALLEEG, EEGTMP, CURRENTSET);' ...
    '  if ~isempty(tmpcom),' ...
    '     EEG = eegh(LASTCOM, EEG);' ...
    '     eegh(tmpcom);' ...
    '     eeglab(''redraw'');' ...
    '  end;' ...
    'end;' ...
    'clear EEGTMP tmpcom;' ];


% call eegplot with the appropriate options
eegplot( EEG.data, 'srate', EEG.srate, 'title', 'Scroll channel activities -- eegplot()', ...
			 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, 'butlabel','LABEL',...
             'events',EEG.event); 

com = [ com sprintf('pop_eeglabel( %s);', inputname(1)) ]; 
return;

