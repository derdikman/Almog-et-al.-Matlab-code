% Rounds x to nearest r
function y = round2r(x,r)
    %y = floor(x) + round((x-floor(x))/r)*r;
    y = round(x/r) * r;