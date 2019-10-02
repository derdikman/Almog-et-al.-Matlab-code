

function [MD_rate_bins, count, time, ang_bins]=Compute_Moving_Directionality...
    (pos_t,pos_x,pos_y,spk_t,pos_vx,pos_vy,parms,ind)
    %this function recieves positional intput and a speed threshold
    %outputs:
    dt=median(diff(pos_t));
    % the moving direction of the animal through out the trail
    pos_MD= atan2(pos_vy,pos_vx);
    spk_MD=interp1(pos_t, unwrap(pos_MD'), spk_t); %get rid of 0, 360 overlap
    spk_MD=wrapToPi(spk_MD);
    n_bins=parms.num_of_direction_bins;
    
    ang_ax =-pi:2*pi/(n_bins-1):pi; %too short
    hist_count_spk = hist(spk_MD,ang_ax);
    hist_time = hist(pos_MD,ang_ax)*dt;

    %NOAM
    ang_ax = -pi:2*pi/(n_bins):pi; %removed bins-1 because ends(pi -pi) are redundant
    count = hist(spk_MD,ang_ax); %spikes in each direction, bin are centers
    count(end) = count(1) + count(end); count(1) = count(end);% noam, because ends are half bins, -180, 180
    time = hist(pos_MD,ang_ax)*dt; %time in each direction
    time(end) = time(1) + time(end); time(1) = time(end);% noam, because ends are half bins, -180, 180
    ang_bins = rad2deg(ang_ax);
    %end NOAM    
    %rate_ang=count./time; %spike rate per direction
    tmp_MD_rate_hist = count./time; %hist_count_spk./hist_time;

    %%convolution with hamming window
%     Win=hamming(10);
%     Win=Win/sum(Win);
%     tmp_MD_rate_hist2=[tmp_MD_rate_hist tmp_MD_rate_hist(2:end) tmp_MD_rate_hist(2:end)]; %nza added(2:end)
%     %tmp_MD_rate_hist2=conv(tmp_MD_rate_hist,Win,'same'); 
%     tmp_MD_rate_hist2=conv(tmp_MD_rate_hist2,Win,'same'); %NOAM
%     % SHAHAF PATCH
%     %this line is meant to solve a problem where the in line 22 the smoothing
%     %causes slightly negative values to apear. 
%     tmp_MD_rate_hist2(tmp_MD_rate_hist2<0)=0;
%     MD_rate_bins=tmp_MD_rate_hist2(length(tmp_MD_rate_hist)+1:2*length(tmp_MD_rate_hist));
    
    MD_rate_bins=tmp_MD_rate_hist;% NOAM IGNORE SMOOTHING
    %MD_rate_hist=tmp_MD_rate_hist((length(Win)/2):length(tmp_MD_rate_hist)-(length(Win)/2));
    
    [max_md_rate,max_bin]=max(MD_rate_bins);

    % find rayleigh score & angle
    norm_val = nansum(MD_rate_bins);
    x_vec = nansum(cos(ang_ax).*MD_rate_bins);
    y_vec = nansum(sin(ang_ax).*MD_rate_bins);
    vec_len = sqrt(x_vec.^2+y_vec.^2);
    MD_rayleigh_score = vec_len/norm_val;
    MD_rayleigh_angle = atan2(y_vec,x_vec); 
    
%     figure;bar(ang_bins,MD_rate_bins);
%     title('MD firing rate per direction deg ');
end
