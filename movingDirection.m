function [rayleigh_angle, phd, shd, rayleigh_score] = movingDirection(c, winms)

dt=median(diff(c.pt)); %winms=10;
st = c.st; pt = c.pt; win = winms;
py1 = c.py(1:end-1); py2 = c.py(2:end);
px1 = c.px(1:end-1); px2 = c.px(2:end);

phd = atan2(py2-py1,px2-px1);
phd=[0; phd];
% the head direction of the animal when spike has ocurred
si = 1; 
if length(pt) > 1 && len(st)>1
   si = discretize(st, [-Inf; mean([pt(1:end-1) pt(2:end)],2); +Inf]);
end
shd = phd(si);

bin_size=deg2rad(1);
radbin = -pi:bin_size:pi;
count = hist(shd,radbin); %spikes per bin
time = hist(phd,radbin)*dt; %time per bin
%consider adding last bins together circular
rate_ang=count./time; %rate per pin

if win > 0
    win=hamming(win);
    win=win/sum(win);
    rate_ang_smoothed=conv([rate_ang rate_ang rate_ang],win,'same');
    rate_ang=rate_ang_smoothed(length(count)+1:2*length(count));
end

x_bin = cos(radbin);
y_bin = sin(radbin);
norm_val = nansum(rate_ang); %denominator to normalize to 1
x_vec = nansum(x_bin.*rate_ang);
y_vec = nansum(y_bin.*rate_ang);
vec_len = sqrt(x_vec.^2+y_vec.^2); %length of x and y component of each bin

rayleigh_score = vec_len/norm_val;
rayleigh_angle=(atan2(y_vec,x_vec));

end