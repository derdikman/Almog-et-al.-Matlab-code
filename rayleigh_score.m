%{ 
https://www.frontiersin.org/articles/10.3389/fphys.2016.00372/full
max mouse running speed 95cm/s
directionality
n=1000; m=20; x=ones(n,1); y=x; xo=1; yo=0;
[a b]=rayleigh_score(1:n,x,y,x+xo,y+yo,x(1:m),y(1:m),x(1:m)+xo,y(1:m)+yo,0);[a rad2deg(b)]
xo=1; yo=0; [1     0]
xo=1; yo=1; [1    45]
xo=0; yo=1; [1    90]
xo=-1;yo=0; [1   180]
xo=0; yo=-1;[1   -90]  
wrapTo360(-90):  270
wrapTo360(-180): 180
wrapTo360(-179): 181
xo=1;yo=sqrt(3): 60
%}
%or (c.pt c.st c.hd) 
function [rayleigh_score, rayleigh_angle ,phd, shd]=rayleigh_score...
    (pt,px1,py1,px2,py2,sx1,sy1,sx2,sy2,winms) %winms originally 10, smoothing hamm window

%rayleigh_score=-1; rayleigh_angle=inf; phd=nan; shd=nan;
%end
%function t
dt=median(diff(pt)); %winms=10;
%head direction
if(nargin == 10) %change to 10? used to be 9
    phd = atan2(py2-py1,px2-px1);
    % the head direction of the animal when spike has ocurred
    shd = atan2(sy2-sy1,sx2-sx1);
    win = winms;
else %rayleigh info has also been calculated, not shd apparently though..
    
    si = 1; st = px1;
    if length(pt) > 1 && len(st)>1
        si = discretize(st, [-Inf; mean([pt(1:end-1) pt(2:end)],2); +Inf]);
    end
    phd=py1;
    shd = phd(si);
    win=10;
    %phd = wrapToPi(pt); old code
    %shd = wrapToPi(px1); old code
    %win = 0;
end


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



%{
   figure;
   subplot(221);bar(radbin,time); title('time spent in each direction');
   subplot(222);bar(radbin,count);title('number of spikes fired in each direction');
   subplot(223);bar(radbin,rate_ang);title('number of spikes fired in each direction');
   subplot(224);title(sprintf('rscore %.3f rangle %.3f',rayleigh_score, rayleigh_angle));
%}

end