%takes RAD as input!
function [mu,s] = meanang(arads,rad)
    arads = toCol(arads);
    mu=circ_mean(arads);
    s=circ_std(arads);
    if nargin ~= 2
        mu =  rad2deg(mu);
        s =  rad2deg(s);
    end
end