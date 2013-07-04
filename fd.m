function [D] = fd(wave, opt)%method, param1,param2)
% fd() - Computed fractal dimension of a waveform
%
% Usage:
%   >> D = fd(X,OPT)
%
% Inputs:
%   X   - (double vector) waveform
%   OPT - (struct) options structure. The following fields are required:
%         .method : 'katz','sevcik','katz_mean','sevcik_mean'
%         .wl     : window length (required for '_mean' methods)
%         .ws     : window shift (required for '_mean' methods)
%
% Outputs:
%   D   - (double) computed fractal dimension
%
%
% References:
%   [1] Katz, M.J., Fractals and the analysis of waveforms, 
%       Comput.Biol.Med. 18: 145, 1988
%   [2] Sevcik, C., A procedure to Estimate the Fractal Dimension of
%       Waveforms, Complexity International, volume 5, 1998, Available online:
%       http://journal-ci.csse.monash.edu.au/ci/vol05/sevcik/
%
%
% Author: German Gomez-Herrero <http://www.cs.tut.fi/~gomezher/index.htm>
%         Institute of Signal Processing
%         Tampere University of Technology, 2009


% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

TOL = 1e-6;

if nargin < 1, help fd; return; end

if ~exist('opt','var'),
    opt = def_fd;
else
    opt = def_fd(opt);
end

method = opt.method;
wl = opt.wl;
ws = opt.ws;

switch(lower(method)),
    case 'katz'
        n = length(wave);
        x = 1:n;
        y = wave;
        % Calculate the diameter
        d = sqrt((x-x(1)).^2+(y-y(1)).^2);
        d = max(d);
        % Calculate the length of the wave
        x = ones(1,(n-1));
        y = wave(2:n)-wave(1:(n-1));
        L = sum(sqrt(x.^2+y.^2));
        D = log10(n)/(log10(d/L)+log10(n));

    case 'sevcik',
        n = length(wave);
        %x = 1:n;
        y = wave;
        % Map the wave to the unit square throught a double linear
        % transformation   
        span = (max(y)-min(y));
        if span < TOL,
            D = 1;
            return; 
        end
        y = (y-max(y))./span;         

        % calculate the length of the wave
        x = (1/(n-1))*ones(1,(n-1));
        y = y(2:n)-y(1:(n-1));
        L = sum(sqrt(x.^2+y.^2));
        D = 1+log(L)/log(2*(n-1));
        
    case 'sevcik_var',
        N = length(wave);        
        %ovlength = wl-ws;
        init = 1:ws:length(wave);
        final = init+wl-1;
        ne = length(init);
        D = zeros(ne,1);
        tmpopt = struct;
        tmpopt.method = 'sevcik';
        for i = 1:ne
           D(i) = fd(wave(init(i):min(final(i),N)),tmpopt); 
        end
        D = var(D);
        
    case 'katz_var',
        N = length(wave);        
        %ovlength = wl-ws;
        init = 1:ws:length(wave);
        final = init+wl-1;
        ne = length(init);
        D = zeros(ne,1);
        tmpopt = struct;
        tmpopt.method = 'katz';
        for i = 1:ne
            D(i) = fd(wave(init(i):min(final(i),N)),tmpopt);
        end
        D = var(D);
        
        
    case 'sevcik_mean',
        N = length(wave);        
        %ovlength = wl-ws;
        init = 1:ws:length(wave);
        final = init+wl-1;
        ne = length(init);
        D = zeros(ne,1);
        tmpopt = struct;
        tmpopt.method = 'sevcik';
        for i = 1:ne
            D(i) = fd(wave(init(i):min(final(i),N)),tmpopt);
        end
        D = mean(D);  
        
     case 'sevcik_window',
        N = length(wave);        
        %ovlength = wl-ws;
        init = 1:ws:length(wave);
        final = init+wl-1;
        ne = length(init);
        D = zeros(ne,1);
        thisopt.method = 'sevcik';
        for i = 1:ne
            D(i) = fd(wave(init(i):min(final(i),N)),thisopt);
        end
        
        
    case 'katz_mean',
        N = length(wave);        
        %ovlength = wl-ws;
        init = 1:ws:length(wave);
        final = init+wl-1;
        ne = length(init);
        D = zeros(ne,1);
        for i = 1:ne
            D(i) = fd(wave(init(i):min(final(i),N)),'katz');
        end
        D = mean(D);


    otherwise
        error('(fd) unknown method %s',method);


end

% subfunction to define the default parameters
% --------------------------------------------
function opt = def_fd(opt)
if nargin < 1 || isempty(opt) || ~isfield(opt,'method'),
    opt.method = 'sevcik';
end
if ~isfield(opt,'wl'),
    opt.wl = [];
end
if ~isfield(opt,'ws'),
    opt.ws = [];
end
