function [pearson_xcov, count_xcor] = compareByMovingDirection(c1, c2, params)
if isfield(params,'lagms')
    nbins = 1;
else
    nbins = params.number_degree_bins; params.lagms = 1000*params.lag_max_secs;
end
if nbins == 1
    pearson_xcov = cell(1,3); count_xcor = cell(1,3);
    [train1,train2] = createMsSpikeTrain(c1.st, 10, c2.st);
    %         spkt1=floor(c1.st*1000)+1; spkt2=floor(c2.st*1000)+1; min_time=min(min(spkt1),min(spkt2)); %TIME IN MS
    %         spkt1 = spkt1-min_time+1; spkt2 = spkt2-min_time+1; max_time=max(max(spkt1),max(spkt2));
    %         train1=zeros(1,max_time);train2 =zeros(1,max_time); train1(spkt1)=1; train2(spkt2)=1; %LOTS OF 0'ss
    
    if params.sigma ~= 0
        %params.sigma = 0.002*1.7^params.sigma; %ORIGINAL
    end
    pearson_xcov{1,1} = timeCorrelationSmoothed(train1,train2,params);
    time_scale=( (1:length(pearson_xcov{1,1})) - ((length(pearson_xcov{1,1})-1)/2) - 1)/1000;
    pearson_xcov{1,2} = time_scale;
    pearson_xcov{1,3} = max(abs(pearson_xcov{1,1}));
    %return
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    px = double(c1.px); py = double(c1.py); pt = double(c1.pt);
    dt = median(diff(pt)); dt = round(dt*1e4)/1e4; % round to 10th of millisec
    vox = zeros(length(pt),1);
    voy = zeros(length(pt),1);
    %figure;subplot(2,2,1); imagesc(c1.rm); title('c1');subplot(2,2,2); imagesc(c2.rm);title('c2');
    %subplot(2,2,3); imagesc(xcorr2(c1.rm, c2.rm)); title('c1xc2'); colormap jet;
    for i = 2:length(vox);
        vox(i)=(px(i)-px(i-1))/(pt(i)-pt(i-1));
        voy(i)=(py(i)-py(i-1))/(pt(i)-pt(i-1));
    end %SMOOTH VELOCITY %test by reconstructing
    [b2,b1] = butter(6,0.06);%0.06);%0.1
    vx = filtfilt(b2,b1,vox);
    vy = filtfilt(b2,b1,voy);
    %REMOVE VELOCITY OVER 3cm/s
    for i = 1:length(vx)
        if sqrt((vx(i)).^2 + (vy(i)).^2) <= 3
            vx(i) = nan; vy(i) = nan;
        end
    end
    %figure();plot(diff(c1.px)/dt);hold on;plot(vx)
    pos_md = atan2(vy,vx);%wrapTo360(rad2deg(atan2d(vy,vx));
    %CHECK IF DELETING A LOT
    c1.st = c1.st(c1.st>=min(pt)&c1.st<=max(pt));
    c2.st = c2.st(c2.st>=min(pt)&c2.st<=max(pt));
    spk1_md = interp1(pt, unwrap(pos_md'),c1.st);
    spk2_md = interp1(pt, unwrap(pos_md'),c2.st);
    %TO 360
    pos_md = wrapTo360(rad2deg(pos_md)); %check??
    spk1_md = wrapTo360(rad2deg(spk1_md)); spk2_md = wrapTo360(rad2deg(spk2_md));
    %BIN BY DIRECTIONS
    bin_edges = (0:360/nbins:360) - 180/nbins; bin_edges(1) = 0; bin_edges(end +1) = 360;% CHECK
    pos_by_bin = discretize(pos_md, bin_edges); pos_by_bin(pos_by_bin==nbins+1) = 1; %conbines two last bins into one, to wrap pi
    %t = abs(diff(pos_by_bin)); %check percent jump in opposite directions bin
    %[sum(t==2) sum(t==0) sum(t==1) sum(t==3) sum(isnan(t))]/length(t)*100;
    %Checked by comparing to pos bins, seems reasonable
    spk1_by_bin = discretize(spk1_md, bin_edges);spk1_by_bin(spk1_by_bin==nbins+1) = 1; %last bin added to first
    spk2_by_bin = discretize(spk2_md, bin_edges);spk2_by_bin(spk2_by_bin==nbins+1) = 1;
    
    %TIME INTERVALS
    time_in_bin = cell(nbins,1);
    i = 1;
    while i <= length(pt)
        cb = pos_by_bin(i);
        increment = true;
        if ~isnan(cb)
            start = pt(i);
            while pos_by_bin(i) == cb && i < length(pt) %will enter first time without changing i
                ends = pt(i);
                i = i+1;
                increment = false;
            end
            start = round(start/dt)*dt; ends = round(ends/dt)*dt; %fancy rounding to dt!
            time_in_bin{cb} = [time_in_bin{cb}; [start:dt:ends]'];
        end
        if increment
            i = i+1;
        end
    end
    
    %BUILD TRAINS
    trains_dir_only=cell(nbins,1);
    for i=1:nbins
        trains_dir_only{i} = zeros(2,length(time_in_bin{i}));
    end
    %merge spikes into time array at closest time within dt
    for i=1:length(spk1_by_bin)
        if ~isnan(spk1_by_bin(i)) %which direction spike in    spike time
            [dif, idx] = min(abs(time_in_bin{spk1_by_bin(i)} - c1.st(i)));
            if dif <= dt
                trains_dir_only{spk1_by_bin(i)}(1, idx) = ... %more than one spike per time unit
                    trains_dir_only{spk1_by_bin(i)}(1, idx) + 1;
            else
                %  fprintf('%d UNSORTED 1 SPIKE %f\n',i,dif);
            end
        end
    end
    for i=1:length(spk2_by_bin)
        if ~isnan(spk2_by_bin(i)) %which direction spike in    spike time
            [dif, idx] = min(abs(time_in_bin{spk2_by_bin(i)} - c2.st(i)));
            if dif <= dt
                trains_dir_only{spk2_by_bin(i)}(2, idx) = trains_dir_only{spk2_by_bin(i)}(2, idx) + 1;
            else
                % fprintf('%d UNSORTED 2 SPIKE %f\n',i,dif);
            end
        end
    end
    binned_trains = trains_dir_only;
    %BIN BY # OF SECS
    bin_size = round(params.time_bin_secs/dt);
    %ALWAYS ZEROS AT END??
    if bin_size > 1
        binned_trains = cell(4,1);
        for i=1:nbins
            z = trains_dir_only{i}';
            if mod(length(z),bin_size) ~= 0
                z = [z;zeros(bin_size-mod(length(z),bin_size),2)];
            end
            %reshapes to rectangle width of bin size, sums cols;
            binned_trains{i}(1,:) = sum(reshape(z(:,1),bin_size, ceil(length(z)/bin_size)));
            binned_trains{i}(2,:) = sum(reshape(z(:,2),bin_size, ceil(length(z)/bin_size)));
        end
    end
    
    %CORRELATION
    pearson_xcov = cell(nbins,3); count_xcor = cell(nbins,3);
    lag_units = round(params.lag_max_secs/params.time_bin_secs);
    for i = 1:nbins
        s1 = binned_trains{i}(1,:);
        s2 = binned_trains{i}(2,:);
        %         % smooth?
        %         win=hamming(5);t1smooth=conv(s1,win,'same');t2smooth=conv(s2,win,'same');
        pearson_xcov{i,1} = xcov(s1,s2,lag_units,'coef')';
        %LOW PASS 0.15%6th order, fc/fs/2 determined empiracally
        %[b2,b1] = butter(6,0.03); pxcsmooth = filtfilt(b2,b1, pearson_xcov{i,1});
        %plot(pxcsmooth); hold on; plot(pearson_xcov{i,1},'y'); legend('smooth','none');
        %pearson_xcov{i,1} = pxcsmooth;
        count_xcor{i,1} = xcorr(s1,s2, lag_units)'; %how many times spike at same point
        %max(max([count_xcor{:,1}]));
        %time scale;
        time_scale=( (1:length(pearson_xcov{i,1})) - ((length(pearson_xcov{i,1})-1)/2) - 1)*params.time_bin_secs;
        pearson_xcov{i,2} = time_scale; count_xcor{i,3}= time_scale;
        %bin edges;
        pearson_xcov{i,3} = [bin_edges(i) bin_edges(i+1)];
        count_xcor{i,3}   = [bin_edges(i) bin_edges(i+1)];
        %plot(pearson_xcov); hold on%figure;plot(count_xcor)
    end
    pearson_xcov{1,3} = [-bin_edges(2) bin_edges(2)]; %correct for 0
    count_xcor{1,3}   = [-bin_edges(2) bin_edges(2)];
    
    fprintf('%.2f ',round(max(max(abs([pearson_xcov{:,1}]))), 2));
    
    
end


%{
            j = 1; k = 0; rem = 1;
        for i = 1:length(spkt1)%params.overlap+1 : length(train1) - paramas.overlap
             t = spkt1(i); %next spk1
             %this spike k in s2 already got deleted and was in overlap of previous s1 spike
             if k > 0 && spkt1(i) < k + params.overlap*1000
                 spkt1(i) = 0;
                 %disp('here');
                 rem = rem+1;%count spks deleted
             end
             %spk2 is within +- overlap         spk2 is smaller than overlap end
             while j <= length(spkt2) && spkt2(j) < t + params.overlap*1000 % in ms
                 %           greater than overlap beginning
                 if spkt2(j) > t - params.overlap*1000
                    %fprintf('%d %d\n',i,j);%spkt1(i),spkt2(j));
                    k = spkt2(j); %last spike in s2 deleted;
                    spkt1(i) = 0; spkt2(j) = 0; %deletes spike in spkt1 and 2 in overlap
                    rem = rem+1;
                 end
                 j = j+1; %next spk2
             end
        end
        spkt1 = spkt1(spkt1~=0); spkt2 = spkt2(spkt2~=0);
%}
