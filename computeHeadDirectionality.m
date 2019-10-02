function [HD_rate_bins, count, time, ang_bins]=computeHeadDirectionality...
    (pos_t,pos_x,pos_y,pos_x2,pos_y2,spk_x,spk_y,spk_x2,spk_y2, bins, ind)

    %Rat head LEDS one on each side <<<
    %dt=pos_t(2)-pos_t(1);
    dt=median(diff(pos_t));
    % the head direction of the animal through out the trail %angle between rat head LEDs <<
    pos_HD = atan2(pos_y2-pos_y,pos_x2-pos_x); %inverse tang gives angle, atan2 avoids ambiguity

    % the head direction of the animal when spike has ocurred
    spk_HD = atan2(spk_y2-spk_y,spk_x2-spk_x);

    n_bins=bins;
    %ang_ax = -pi:2*pi/(n_bins-1):pi;
    ang_ax = -pi:2*pi/(n_bins):pi; %removed bins-1 because ends(pi -pi) are redundant
    count = hist(spk_HD,ang_ax); %spikes in each direction, bin are centers %binning
    count(end) = count(1) + count(end); count(1) = count(end);% noam, because ends are half bins, -180, 180
    time = hist(pos_HD,ang_ax)*dt; %time in each direction
    time(end) = time(1) + time(end); time(1) = time(end);% noam, because ends are half bins, -180, 180
    HD_rate_bins=count./time; %spike rate per direction
    
    ang_bins = rad2deg(ang_ax);
    
    %rayleigh_score
    x_val = cos(ang_ax);
    y_val = sin(ang_ax);
    norm_val = nansum(HD_rate_bins);
    x_vec = nansum(cos(ang_ax).*HD_rate_bins);
    y_vec = nansum(sin(ang_ax).*HD_rate_bins);
    vec_len = sqrt(x_vec.^2+y_vec.^2);
    rayleigh_score = vec_len/norm_val;
    rayleigh_angle=(atan2(y_vec,x_vec));

    close all;
%     figure;
%     subplot(1,3,1);bar(ang_bins,HD_rate_bins); title('HD firing rate per direction deg');
%     subplot(1,3,2);bar(ang_bins,time);       title('time spent in each direction');
%     subplot(1,3,3);bar(ang_bins,count);      title('number of spikes fired in each direction');
%     figure;bar(ang_ax,rate_ang); title('firing rate in each direction');%direction');
    %t = [22,23,45,67,68, 158,190, -158], hist(t,aa)
end

