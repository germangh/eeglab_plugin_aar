function [flx] = fl(x,fp_t)

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com

fp_U=30;
fp_L=-30;
fp_method = 'R';
if fp_t==Inf, flx=x; return; end
flx = zeros(size(x));

index = find(x~=0);
x = x(index);

for i = 1:length(index),
    exponent = 0;
    mantissa = abs( x(i) );
    while mantissa >= 1
        mantissa = mantissa/10;
        exponent = exponent+1;
    end
    while mantissa < .1
        mantissa = mantissa * 10;
        exponent = exponent-1;
    end

    % Round or truncate, as indicated by fp_method
    if fp_method == 'R'
        mantissa = round( mantissa*(10^fp_t) ) / (10^fp_t);
    elseif fp_method == 'T'
        mantissa = floor( mantissa*(10^fp_t) ) / (10^fp_t);
    end

    if mantissa == 1
        mantissa = mantissa/10;
        exponent = exponent+1;
    end

    if exponent > fp_U
        disp( 'WARNING overflow' )
    elseif exponent < fp_L
        disp( 'WARNING underflow' )
    end

    flx(index(i)) = sign(x(i)) * mantissa * 10^exponent;
end