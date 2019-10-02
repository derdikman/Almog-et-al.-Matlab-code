%%%%%%%%%
% N Z A %
% grid cells: 228
% pairs: 749
%%%%%%%%%

%{
** TO DO

%}

function musciomol()
    warning off;
    dl = 'C:\Noam\Data\muscimol\DB_musc_MEC\';
    ds = 'C:\Noam\Data\muscimol\noam_output\';
    %%%
    binsize = inf;  %change to 45   WAS 45
    fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh.mat',binsize);
    %fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_MINI.mat',binsize);
    %    tic %1800 sec   %LOAD FROM FILES
    %%{
    tic;
    files = dir(strcat(dl,'DB*.mat'));
    for i = 1:length(files) %CHECK 36
        %i = 81;
        data{i} = load(strcat(dl,files(i).name),'db');
        data{i}.db.ind = i;
        cells{i} = process(data{i}.db, binsize); 
        fprintf('%d %.1f \n',i, (i*100/length(files)));
    end
    toc
    
    tic
    save(fn,'cells');
    toc
    %}
    
    
    fprintf('loading %s ',fn); %ascii 48
    tic; cells = load(fn); toc;
    %c1=cells{71}.before;, c2=cells{75}.before; n = 10;
    %batchTimeCorrelations (c1, c2, n);
    %stop
    
    %
    
    
    
    t = [];
    for i = 1:length(cells)
        c = cells{i}.before;
        t(end+1) = c.max_r;
    end
    [m,i] = max(t);
    c = cells{i}.before;
    max(diff(c.pt))
    STOP
    
    
    %saveSpikesByDirection(params, 8);
    
    %driftWindowMain();
    %stop
    
    
    %plotDirectionMain();
    

    %stop
    [groups cells] = findSimultaneouslyRecordedCells(cells);
    %k = fieldnames(groups);
   
    disp('done and done'); 
    stop
    %TEMPORAL CORRELATION
    %{
    figure();
    tic
    for i = 1:length(k)
        %disp(i);
        colormap jet;
        g = groups.(k{i});
        [mc, lc] = plotGroupTemporralCorrs(g, 500);
%         subplot(10,9,2*i-1);
%         imagesc(lc./max(lc(:))); axis off; axis equal; caxis([-1 1]);
%         title(sprintf('i%d:[%.1f %.1f]',i, min(lc(:)), max(lc(:))));
%         set(gca,'FontSize',8)
%         subplot(10,9,2*i);
%         imagesc(mc); axis off; axis equal; %caxis([-1 1]);
%         [m, ind] = max(mc(:));
%         title(sprintf('i%d:[%.1f %.3f]',i, lc(ind), m));
%         set(gca,'FontSize',8)
         subplot(5,8,i);
         imagesc(mc,[0,0.5]); axis off; axis equal; %caxis([-1 1]);
         [m, ind] = max(mc(:));
         title(sprintf('i%d(%.3f)(%.1f)',i, m, lc(ind)));
         set(gca,'FontSize',9); 
        
        fprintf('%d %.1f %.3f\n',i,lc(ind),m);
    end
    set(suptitle('Highest Temporal-Correlation/Lag for Each Pair in Group, Normalized'),'Interpreter', 'none');
    toc; disp('time corr');
    %}
    
    
    %i = 5; plotGroup(groups.(k{i}), k{i}, i, binsize, midmax);
    
    %plot groups
    %args.cells = cells;
    args.binsize = 10;
    args.groups = groups;
    args.grid_thresh = 0.5;
    plotGroup(args, 1);
    stop;


    %print group properties
    for i = 5%1:length(k);
        g = groups.(k{i});
        maxRates = g(1); maxRates = length(maxRates.middle); %t = length(g);
        %t.a; length(t.after.px) .exists length(t.middle)
        fprintf('%s * %d \n',k{i} ,maxRates);
        %fprintf('%s %d', k{i}, length(g));
        for j = 1:length(g) %in group
            tt = g(j); tt = length(tt.middle);
            %.a length(tt.after.px) .exists length(tt.middle)
            if maxRates ~= tt               %not(strcmp(t,tt))
                fprintf('%d ', tt);
                maxRates = tt;
            end; %fprintf('%d %d\n', g(j).tet, g(j).cel);
        end; disp('*'); 
    end



    %pairs
    tic
    pairs = find_pairs(cells);%find_pairs_of_cells_mosimol(data);
    d = 47; %65 23 86 235
    for i = 48:49 %47:length(pairs)
        r1i = pairs(i,1); r2i = pairs(i,2);
        fprintf('%.1f of %d ',i*100/length(pairs),length(pairs));
        r1 = cells{r1i};%prosess(data{r1i}.db);
        r2 = cells{r2i};%prosess(data{r2i}.db);
        %plotPair(r1, r2, i);
    end
    toc
    disp('muscimol done');
    %{
        - Gilad's paper
        - to what extent do cells disappear together
        - max simul recorded
        - module
    %}
end



function driftWindowMain()
    fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_a.mat', 10);
    fprintf('loading **bin size**%d*\n',10); 
    tic; cells = load(fn); cells = cells.cells; toc;
    groups = findSimultaneouslyRecordedCells(cells);
    k = fieldnames(groups);   
    for i = 1:length(k)
        g= groups.(k{i});
        good = [g(1)]; bad = [g(1)]; 
        for j = 1:length(g) %put low grid scores at end
            if g(j).before.gridscore > 0.5; %GRID THRESH
                good(end+1) = g(j);
            else
                bad(end+1) = g(j);
            end; 
        end
        bad = bad(2:end); good = good(2:end);
        %iterate through group
        for j = 1:length(good)
            c = good(j);
            s = c.before;
            driftWindow(s,60);
            
     set(g,'PaperPositionMode', 'auto');
    print(g,sprintf('temp/%s.png',titl),'-dpng', '-r0');
    close(g);           
            
    m.rm = m.rm ./ length(windowed);
    g = figure(); colormap jet;
    subplot(121);
    imagesc(m.rm); set(gca,'ydir','normal'); axis square;
    subplot(122);
    m.ac = Cross_Correlation(m.rm, m.rm);
    imagesc(m.ac); set(gca,'ydir','normal'); axis square;
    titl = sprintf('mean rate map win %ds.png',winsecs);
    suptitle(titl);
    set(g,'PaperPositionMode', 'auto');
    print(g,sprintf('temp/%s.png',titl),'-dpng', '-r0');
    close(g);

%set(h, 'PaperPosition', [0 0 300 100*r]);
set(h, 'paperunits', 'inches');
set(h, 'PaperSize', [4 r*10]);
set(h,'PaperPositionMode', 'manual'); set(h,'PaperSizeMode', 'manual'); 
set(h,'PaperSizeMode', 'auto');
print(g,sprintf('temp/%stest3.png',''),'-dpng', '-r0');    
stop
print('test.pdf');
export_fig test2.pdf
print('test.png');
            
        end
    end
    
    
    
end

function plotDirectionMain()
    params.bindeg = 15;
    params.lag  = 2.4;
    params.binspike  = 0.06;
    params.sigma = 2;
    params.version = 'n';
    params.grid_thresh = 0.5;
    params.sesh = 'before';
    fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_a.mat', 10);
    fprintf('loading **bin size**%d*\n',10); %ascii 48
    tic; cells = load(fn); cells = cells.cells; toc;
    params.groups = findSimultaneouslyRecordedCells(cells);
    k = fieldnames(params.groups);
    for i = 1:length(k)
        params.gid = i;
        [h,v] = plotByDirectionMain(params);
        if ~isempty(v)
            sdir = 'C:\Noam\Output\muscimol\HDMD\'; debug = '';
            filename = sprintf('%s%s%s.png',sdir, debug, v); disp(filename);
            %stop above before save
            h.PaperPositionMode = 'auto'; print(filename, '-dpng','-r0');
            close(h);
        end
    end
    stop
end

function saveSpikesByDirection(params, degbins)
    tic; files = dir(strcat(params.dir_load,'DB*.mat'));
    prnt = true;
    for i= 1:length(files)
        fprintf('%d %.1f \n',i, (i*100/length(files)));
        a = load(strcat(params.dir_load,files(i).name));
        %a = load(strcat(params.dir_load,'DB_MUSC_MEC_001_11468_010306_t2c1.mat '));
        p = a.db.B(1).pos_data; s = a.db.B(1).spike_data;

        %HEAD DIRECTION &&&&&&
        % (pos_t,pos_x,pos_y,pos_x2,pos_y2,spk_x,spk_y,spk_x2,spk_y2, bins)
        [HD_rate_bins,HD_count, HD_time, ang_bins]= computeHeadDirectionality...
            (p.t,p.x1,p.y1,p.x2,p.y2,s.x,s.y,s.x2,s.y2, degbins);

        %MOVING DIRECTION
        assert(length(s.x) == length(s.ts), 'length(s.x) == length(s.ts)');
        a = rat_trajectory(p.x1, p.x2, p.y1, p.y2, p.t, s.x, s.x2, s.y, s.y2, s.ts);
        %tic; cells = load(fn); cells = cells.cells; toc; a = cells{1}.before;
        px = double(a.px); py = double(a.py); pt = a.pt; st = a.st;
        % computing velocity
        parms.smoothing_parameter_for_velocity = 10e-3; %10e-5;
        tx = csaps( 1:length(pt),px);%parms.smoothing_parameter_for_velocity );
        ty = csaps( 1:length(pt),py);%parms.smoothing_parameter_for_velocity ); %param = 10e-5 or 10e-3
        %spline and smoothing, off by one binning and smoothing
        % Compute the velocity via a derivative of the pp-form of the spline (fnder.m):
        dt=median(diff(pt));
        vx = fnval( fnder(tx),1:length(pt))/dt; %tx
        vy = fnval( fnder(ty),1:length(pt))/dt; %ty
        %vx = diff(px)/dt; vx = [0 vx];
        %vy = diff(py)/dt; vy = [0 vy];
        bins = 8; parms.num_of_direction_bins = bins;
        [MD_rate_bins,MD_count,MD_time, ang_bins]=...
            computeMovingDirectionality(pt,px,py,st,vx,vy,parms); %Compute_Moving_Directionality
        
        %correct
        phd = atan2(p.y2-p.y1,p.x2-p.x1);
        pmd = atan2(vy,vx)';
        cellByDeg{i}.ind = i; cellByDeg{i}.ang_bins = ang_bins;
        cellByDeg{i}.HD_rate = HD_rate_bins;  cellByDeg{i}.MD_rate = MD_rate_bins;
        if prnt
            h = figure('Position', [0, 0, 1000, 500]);%(2*m+5)*120 , 140*r + height]); %
            set(gca,'LooseInset', get(gca,'TightInset'));colormap jet;
            titl = sprintf('Head_Moving_Direction_%dbins_i%d', bins, i);
            subplot(2,3,1);bar(ang_bins,HD_rate_bins); title('HD firing rate');
            subplot(2,3,2);bar(ang_bins,HD_time);       title('Time');
            subplot(2,3,3);bar(ang_bins,HD_count);      title('Spikes');     
            subplot(2,3,4);bar(ang_bins,MD_rate_bins); title('MD firing rate');
            subplot(2,3,5);bar(ang_bins,MD_time);       title('Time');
            subplot(2,3,6);bar(ang_bins,MD_count);      title('Spikes');        
            ax = findobj(gcf,'Type','Axes');
            for q=1:length(ax)
                set(ax(q),'FontSize',9);%axis(ax(q),'equal');%axis(ax(q),'off');%title(ax(i),{'Very Nice'})
            end
            set(suptitle(titl),'Interpreter', 'none');
            sdir = 'C:\Noam\Output\muscimol\HDMD\'; debug = '';
            filename = sprintf('%s%s%s.png',sdir, debug, titl); disp(filename);
            h.PaperPositionMode = 'auto'; print(filename, '-dpng','-r0');
            close(h)
        end
    end
    %SAVE VAR
    fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\direction_spikes_%ddegs.mat',degbins);
    disp(fn);
    save(fn,'cellByDeg');
end



%{
    % correlation bt 2 simultaneously recorded cells, same module(?) not overlapping(?)
    % for each 8 directions, temporal correlation between the 2 cells, vary amount of time
    % use moving direction
    % correction not for this
    % correction cell by cell
 - Take group
 - Find grid cells
 - Take one cell, for each direction, try to correlate each other cell
 - for each group x for each cell x for each direction x all other cells
 - for g x c x d x cs
 - COMPARE: 75 76, 226 228 
%}


%compares moving bins of c1 to entire c2




function correlation_sanity_check(c1, c2)
    nbins = 4;
    %MOVING DIRECTION
    smooth = 10e-7;%10e-3 %10e-5; %Empiracally so bins don't change rapidly for 8 bins!
    tx = csaps( 1:length(c1.pt),double(c1.px), smooth);%10e-3; %10e-5;
    ty = csaps( 1:length(c1.pt),double(c1.py), smooth);
    dt = median(diff(c1.pt));
    vx = fnval( fnder(tx),1:length(c1.pt))/dt; 
    vy = fnval( fnder(ty),1:length(c1.pt))/dt;
    pos_md = atan2d(vy,vx);%wrapTo360(rad2deg(atan2d(vy,vx)));    
    %TO 360 
    pos_md = wrapTo360(rad2deg(pos_md)); %check??
    %BIN BY DIRECTIONS
    bin_edges = (0:360/nbins:360);% - 180/nbins; bin_edges(1) = 0; bin_edges(end +1) = 360; 
    pos_by_bin = discretize(pos_md, bin_edges)';% pos_by_bin(pos_by_bin==nbins+1) = 1; %conbines two last bins into one, to wrap pi
    a = diff(c1.pt);
    unique(a)
    edges = round(edges);
    direction_array = interp1(c1.pt, unwrap(pos_md'),edges);   
    %SPIKES
    c1.st = c1.st(c1.st>=min(c1.pt)&c1.st<=max(c1.pt)); %removes spikes out of bound pos times
    c2.st = c2.st(c2.st>=min(c1.pt)&c2.st<=max(c1.pt)); %removes spikes out of bound pos times
    %to 360
    spk1_md = interp1(c1.pt, unwrap(pos_md'),c1.st);
    spk2_md = interp1(c1.pt, unwrap(pos_md'),c2.st);
    spk1_md = wrapTo360(rad2deg(spk1_md)); spk2_md = wrapTo360(rad2deg(spk2_md));
    %FIRING RATE
    spkt1 = c1.st; spkt2 = c2.st;
    spkt1 = floor(spkt1*1000)+1; spkt2=floor(spkt2*1000)+1; %in MILLISECS
    max_time=max(max(spkt1),max(spkt2));
    %firing rate;
    bin_msecs = 100;
    firing1 = zeros(max_time,1); firing2 = zeros(max_time,1);
    firing1(spkt1) = spkt1;  firing2(spkt2) = spkt2;    
    [firing1, edges] = histcounts(firing1,round(max_time/bin_msecs));%(x,nbins) 
    firing2 = histcounts(firing2,round(max_time/bin_msecs));%now in units of binsecs
    firing1(1) = 0; firing2(1) = 0; %spikes in first bin_msecs ignored
    
    %   
    spk1_by_bin = discretize(spk1_md, bin_edges);  %spk1_by_bin(spk1_by_bin==nbins+1) = 1;
    spk2_by_bin = discretize(spk2_md, bin_edges);    
    
    [P,F] = pwelch(y);
    helperFilterIntroductionPlot1(F,P,[60 60],[-9.365 -9.365],...
  {'Original signal power spectrum', '60 Hz Tone'})
    fc = 500;
    fs = 10000;
    figure;
    plot(y); hold on;
    plot(filter(b,a,t));
    [z,p,k] = butter(10,0.05);
    sos = zp2sos(z,p,k);
    grpdelay(sos,128);
    
    
end

%spkt1 = g(i).before.st;
function [pxcsmooth, count_xcor] = time_correlation(spkt1,spkt2, params)
    %               generating trains
    %scales times x 1000 (ms), +1 so non-zero
    bin_secs = params.time_bin_secs;
    lag = params.lag_max_secs/bin_secs;
    spkt1=floor(spkt1*1000)+1; spkt2=floor(spkt2*1000)+1; %in MILLISECS
    min_time=min(min(spkt1),min(spkt2));
    spkt1 = spkt1-min_time+1; spkt2 = spkt2-min_time+1;%normalizes
    max_time=max(max(spkt1),max(spkt2));
    %generating the spike train
    train1=zeros(1,max_time);train2 =zeros(1,max_time); 
    train1(spkt1)=1; train2(spkt2)=1;%array at these indices(time*1000) will be 1, else 0
    %[train1,train2,count_xcor,pearson_xcor] = My_Xcor(train(1,:),train(2,:));
    assert(lag>=1,'lag=>1');
    bin_msecs = bin_secs *1000;
    %binned!
    train1 = histcounts([1 spkt1 max_time],round(max_time/bin_msecs));%(x,nbins) 
    train2 = histcounts([1 spkt2 max_time],round(max_time/bin_msecs));%now in unites of binsecs
    %SUBTRACT FIRST AND LAST ADDED SPIKE
    train1(1) = train1(1)-1;train2(1) = train2(1)-1;train1(end) = train1(end)-1;train2(end) = train2(end)-1;
    win=hamming(5); %light 3, could make 5)
    train1smooth=conv(train1,win,'same');
    train2smooth=conv(train2,win,'same');
    %generating the cross corelations normelized to pearson
    pearson_xcor=xcov(train1smooth,train2smooth,lag,'coef')'; %500 parms.max_lag NEED THIS <<<
	[b,a] = butter(6,0.03); %LOW PASS 0.15%6th order, fc/fs/2 determined empiracally 
    pxcsmooth = filtfilt(b,a,pearson_xcor);
    %close all; plot(pearson_xcor); hold on; plot(pxcsmooth);plot(xcov(train1,train2,lag,'coef')'); legend('smooth','none');
    count_xcor=xcorr(train1smooth,train2smooth, lag)'; %how many times spike at same point <<<
end


function [maxcorr, lagcorr] = plotGroupTemporralCorrs(g,lag)
    maxcorr = zeros(length(g));
    lagcorr = zeros(length(g));
    for i=1:length(g) 
        t1 = g(i).before.st;
        for j=i+1:length(g)
            t2 = g(j).before.st;
            [p, c] = time_correlation(t1,t2);
            %A(i,j) = corr(t1, t2);           
            [m, ind] = max(p);
            maxcorr(i,j) = m;
            lagcorr(i,j) = ind - (lag + 1); %+1 so that 1 - (500 + 1) = -500
        end
    end  
end


%{
remove corner cases
%}



function plotPair(r1, r2, id)
    figure('Position', [100, 100, 1200, 900]);
    set(gca,'LooseInset',get(gca,'TightInset'));
    colormap jet; bad = '';
    rmt = 0.1; % Rate Map Threshold
    %%% set(gca,'YDir','normal') instead of flip???
    
    %n = min(length(r2.middle), length(r1.middle));
    if length(r2.middle) ~= length(r1.middle)
        %fprintf('%d: error plotting, middle sessions different lengths', id);%bad = '*';
    end

    set(suptitle(sprintf('%sP%d(i%d,i%d): rat %s %s C1(t%dc%d %s %s) C2(t%dc%d %s %s)',...
        bad,id,r1.ind,r2.ind,r1.id,r1.date,r1.tet,r1.cel,r1.a,r1.type,r2.tet,r2.cel,r2.a,r2.type)),...
        'Interpreter', 'none');

    %line up bin times
    i = 1; i1 = 1; i2 = 1; mid1 = {}; mid2 = {};
    while i1 <= length(r1.middle) && i2 <= length(r2.middle)
        if r1.middle{i1}.pt(1) == r2.middle{i2}.pt(1) %Will only show traj,ratemap/ac for good ratemaps
            if not(r1.middle{i1}.max_r < rmt || r1.middle{i1}.max_r == 50 ...
                    || r2.middle{i2}.max_r < rmt || r2.middle{i2}.max_r == 50)
                mid1{i} = r1.middle{i1};
                mid2{i} = r2.middle{i2};
                i = i + 1;
            end
            i1 = i1 + 1;
            i2 = i2 + 1;
        elseif r1.middle{i1}.pt(1) < r2.middle{i2}.pt(1)
            i1 = i1 + 1;
        else
            i2 = i2 + 1;
        end
    end
    r1.middle = mid1;
    r2.middle = mid2;
    m = length(mid1);
    n = m +1+1; %how many cols in figure, +1 before +1 mcross
    after = 0; %display after;
    if not(r1.after.max_r < rmt || r1.after.max_r == 50 ...
            ||r2.after.max_r < rmt || r2.after.max_r == 50)
        n = n + 1;
        after = 1;
    end
    r = 8; %Rows for figure (trajectory ratemap autoc xcross axcross)


    %%%plot trajectory
    subplot(r, n, 1 + n*0);
    plot(r1.before.px, flip(r1.before.py));
    hold on
    scatter(r1.before.sx, flip(r1.before.sy), '.'),...
        xlim([0 100]), ylim([0 100]), ...
        title(sprintf('c1: before %d',length(r1.before.st)));
    axis off; axis equal;
    subplot(r, n, 1 + n*1);
    plot(r2.before.px, flip(r2.before.py));
    hold on
    scatter(r2.before.sx, flip(r2.before.sy), '.'),...
        xlim([0 100]), ylim([0 100]), axis off; axis equal;
    title(sprintf('c2: %d',length(r2.before.st))); %TITLE
    for i = 1:m;
        %if length(r1.middle{i}.px) > 1 && length(r2.middle{i}.px) > 1 %Display if more than one spike
        subplot(r, n, i+1 + n*0);
        plot(r1.middle{i}.px, flip(r1.middle{i}.py));
        hold on
        scatter(r1.middle{i}.sx, flip(r1.middle{i}.sy), '.'); xlim([0 100]), ylim([0 100]);
        %title(sprintf('%.0f-%.0f(%d)', min(r1.middle{i}.pt)/60,...
        %max(r1.middle{i}.pt)/60, length(r1.middle{i}.st)));
        title(sprintf('%.0f-%.0f(%d)', r1.middle{i}.pt(1)/60,...
            r1.middle{i}.pt(length(r1.middle{i}.pt))/60, length(r1.middle{i}.st)));
        axis off; axis equal;
        subplot(r, n, (i+1) + n*1);
        plot(r2.middle{i}.px, flip(r2.middle{i}.py));
        hold on
        scatter(r2.middle{i}.sx, flip(r2.middle{i}.sy), '.'), xlim([0 100]), ylim([0 100]);
        %title(sprintf('%d',length(r2.middle{i}.st)));
        title(sprintf('%.0f-%.0f(%d)', r2.middle{i}.pt(1)/60,...
            r2.middle{i}.pt(length(r2.middle{i}.pt))/60, length(r2.middle{i}.st)));
        axis off; axis equal;
    end
    %if length(r1.after.px) > 1 && length(r2.after.px) > 1
    if after
        subplot(r, n, n*1-1);
        plot(r1.after.px, flip(r1.after.py));
        hold on
        scatter(r1.after.sx, flip(r1.after.sy), '.'),...
            xlim([0 100]), ylim([0 100]), axis off; axis equal;
        title(sprintf('after %d',length(r1.after.st)));
        subplot(r, n, n*2-1);
        plot(r2.after.px, flip(r2.after.py));
        hold on
        scatter(r2.after.sx, flip(r2.after.sy), '.'),...
            xlim([0 100]), ylim([0 100]); axis off; axis equal;
        title(sprintf('%d',length(r2.after.st)));
    end

    %%%plot rm
    subplot(r, n, 1 + n*2);
    imagesc(r1.before.rm), title(sprintf('c1: %.0fHz', r1.before.max_r)); axis off; axis equal;
    subplot(r, n, 1 + n*3);
    imagesc(r2.before.rm), title(sprintf('c2: %.0fHz', r2.before.max_r)); axis off; axis equal;
    for i = 1: m;
        subplot(r, n, i+1 + n*2);
        imagesc(r1.middle{i}.rm);title(sprintf('%.0fHz', r1.middle{i}.max_r));axis off; axis equal;
        subplot(r, n, i+1 + n*3);
        imagesc(r2.middle{i}.rm);title(sprintf('%.0fHz', r2.middle{i}.max_r));axis off; axis equal;
    end                           %after
    if after
        subplot(r, n, n*3-1);
        imagesc(r1.after.rm);title(sprintf('%.0fHz', r1.after.max_r));axis off; axis equal;
        subplot(r, n, n*4-1);
        imagesc(r2.after.rm);title(sprintf('%.0fHz', r2.after.max_r));axis off; axis equal;
    end
    % mcross
    t1{1} = r1.before.rm; t2{1} = r2.before.rm;
    for i = 1:length(mid1)
        t1{i+1} = mid1{i}.rm; t2{i+1} = mid2{i}.rm;
    end
    if after
        t1{m+2} = r1.after.rm; t2{m+2} = r2.after.rm;
    end
    subplot(r, n, n*3);
    imagesc(allCorr(t1,t1));title('cross per');axis off; axis equal;
    subplot(r, n, n*4);
    imagesc(allCorr(t2,t2));title('cross per');axis off; axis equal;

    %%%plot ac
    subplot(r, n, 1 + n*4);
    imagesc(r1.before.ac), title(sprintf('c1: %.2f', r1.before.gridscore));
    axis off; axis equal;
    subplot(r, n, 1 + n*5);
    imagesc(r2.before.ac), title(sprintf('c2: %.2f', r2.before.gridscore));
    axis off; axis equal;
    for i = 1: m;
        subplot(r, n, i+1 + n*4);
        imagesc(r1.middle{i}.ac); title(sprintf('%.2f', r1.middle{i}.gridscore)); axis off; axis equal;
        subplot(r, n, i+1 + n*5);
        imagesc(r2.middle{i}.ac); title(sprintf('%.2f', r2.middle{i}.gridscore)); axis off; axis equal;
    end
    if after
        subplot(r, n, n*5-1);
        imagesc(r1.after.ac); axis off; axis equal;
        title(sprintf('%.2f', r1.after.gridscore));
        subplot(r, n, n*6-1);
        imagesc(r2.after.ac); axis off; axis equal;
        title(sprintf('%.2f', r2.after.gridscore));
    end
    % mcross
    t1 = {}; t2 = {}; t1{1} = r1.before.ac; t2{1} = r2.before.ac;
    for i = 1:length(mid1)
        t1{i+1} = mid1{i}.ac; t2{i+1} = mid2{i}.ac;
    end
    if after
        t1{m+2} = r1.after.ac; t2{m+2} = r2.after.ac;
    end
    subplot(r, n, n*5);
    imagesc(allCorr(t1,t1)); title('cross per');axis off; axis equal;
    subplot(r, n, n*6);
    imagesc(allCorr(t2,t2));title('cross per');axis off; axis equal;

    %%% cross corr
    subplot(r, n, 1 + n*6);
    cc = Cross_Correlation(r1.before.rm, r2.before.rm);
    imagesc(cc); title(sprintf('Xcorr: %s',''));axis off; axis equal;
    %acc
    subplot(r, n, 1 + n*7);
    acc = Cross_Correlation(cc, cc); acc = getSubmatrixFromCenter(acc, cc);
    imagesc(acc);
    title(sprintf('AcXcorr: %.2f', gridscore(acc, -1))); axis off; axis equal;
    for i = 1: m;
        subplot(r, n, i+1 + n*6);
        cc = Cross_Correlation(r1.middle{i}.rm, r2.middle{i}.rm);
        imagesc(cc);axis off; axis equal;
        %acc
        subplot(r, n, i+1 + n*7);
        acc = Cross_Correlation(cc, cc); acc = getSubmatrixFromCenter(acc, cc);
        imagesc(acc);axis off; axis equal;title(sprintf('%.2f', gridscore(acc, -1)));
    end

    if after
        subplot(r, n, n*7-1);
        cc = Cross_Correlation(r1.after.rm, r2.after.rm); imagesc(cc);axis off; axis equal;
        %acc
        subplot(r, n, n*8-1); acc = Cross_Correlation(cc, cc);
        getSubmatrixFromCenter(acc, cc); imagesc(acc);
        title(sprintf('%.2f', gridscore(acc, -1)));axis off; axis equal;
    end




    sdir = 'C:\Noam\Output\muscimol\pairs\';
    filename = sprintf('%s%d_Rat_%s_date_%s_C%d_t%d_c%d_C%d_t%d_c%d.png',sdir, id, r1.id, r1.date, r1.ind, r1.tet, r1.cel, r2.ind, r2.tet, r2.cel);
    disp(filename);
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(filename, '-dpng','-r0');
    close;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end plot
function f = flip(x)
    m = max(x);
    f = m-x;
end

%not perfect, submatrix will be a dimension higher than B, if odd
%(could fix using length instead of round(B)
function S = getSubmatrixFromCenter(A, B) %A big matrix, B little matrix
raca = round(size(A)/2);
rbcb = round(size(B)/2);
S = A((raca(1)-rbcb(1)+1:raca(1)+rbcb(1)),((raca(2)-rbcb(2)+1:raca(2)+rbcb(2))));
end

%refactored
% function s = makeSession(p,s,n,mnxt)
%     s = rat_trajectory(p.x1, p.x2, p.y1, p.y2, p.t, s.ts, mnxt); %s.x, s.x2, s.y, s.y2,
%     [s.rm, s.max_r] = Create_Rate_Map(s.px, s.py, s.pt, s.sx, s.sy, s.st, true, n);
%     %s.acOrig = Cross_Correlation(s.rm, s.rm); 
%     s.ac = xcorr2(s.rm); 
%     s.gridscore = gridscore2(s.ac, 2);
%     s.module = Find_Module(imgaussfilt(s.ac, 3,'FilterDomain','spatial'));
%     s.exists = false;
%     if s.max_r ~= 0 %&& s.gridscore ~= -2 %&& s.max_r ~= 50 
%         s.exists = true;
%     end
% end



%refactored
% function c = rat_trajectory(px1, px2, py1, py2, pt, st, mnxt) %sx1, sx2, sy1, sy2,
%     pi = mnxt(1) <= pt & pt <= mnxt(2); si = mnxt(1) <= st & st<=mnxt(2);
%     pt = double(toCol(pt(pi))); st = double(toCol(st(si)));
%     minx = min(min(px1),min(px2));miny = min(min(py1),min(py2));
%     px1 = double(toCol(px1(pi)) - minx +.00001); 
%     py1 = double(toCol(py1(pi)) - minx +.00001);
%     px2 = double(toCol(px2(pi)) - miny +.00001); 
%     py2 = double(toCol(py2(pi)) - miny +.00001);
%     
%     %PATCH TRAJECTORY
%     dt = median(diff(pt));
%     p1 = patchTrajectoryLinear(pt,px1,py1,dt,1.1*dt); 
%     p2 = patchTrajectoryLinear(pt,px2,py2,dt,1.1*dt); 
%     c.pt = p1.t; %PT%
%     
%     %case for no spikes
%     c.st = st;
%     if isempty(st)
%         st = c.pt(1);
%     end
%     
%     %get closest time in px
%     si = 1; 
%     if length(c.pt) > 1
%         si = discretize(st, [-Inf; mean([c.pt(1:end-1) c.pt(2:end)],2); +Inf]);
%     end
%     
%     %base sx sy on closest time of spike to data in px py
%     [c.rayleigh_score, c.rayleigh_angle, c.hd] =...
%         rayleigh_score(c.pt,p1.x,p1.y,p2.x, p2.y, p1.x(si),p1.y(si),p2.x(si),p2.y(si));
%     
%     c.px =  mean([p1.x, p2.x], 2);
%     c.px = c.px - min(c.px) + 0.00001; %no zeros
%     c.py =  mean([p1.y, p2.y], 2);
%     c.py = c.py - min(c.py) + 0.00001;
%     c.st = st;
%     %base sx sy on closest time of spike to data in px py
%     c.sx = c.px(si);
%     c.sy = c.py(si);
%     
% %     c.sx = mean([sx1, sx2], 2); c.sx = c.sx - min(c.sx) + 0.00001;
% %     c.sy = mean([sy2, sy2], 2); c.sy = c.sy - min(c.sy) + 0.00001;
%     %{figure();plot(c.px, c.py);hold on; plot(c.sx, c.sy,'.');%}
% end

function test()
    sa=0, sb= 0;
    for i = 1:length(bins)
        fprintf('%d: min %.3f max %.3f \n',i, min(bins{i}.pt)/60, max(bins{i}.pt)/60);
        sa = sa + length(bins{i}.pt);
        sb = sb + length(bins{i}.st);
    end
    fprintf('%d  %d \n',sa, sb);
end

%tr = trajectory; bs = bin size in seconds; mt = max time
function bins = bin_trajactory(tr, bs, mt)%
minSpikes = 1;%min  spikes to be a bin
minBinLength = 1; %in minutes
bins{1,1} = []; b = 1;
ind = 0; %index in current bin
t = tr.pt; st = tr.st;
%t0 = ceil(tr.pt(1)/(5*60))*5*60 - bs; %statrs at 5min intervals %tr.pt(1);
t0 = tr.pt(1);
%SHOULD REWRITE VECTORIZED
for i = 1: length(t) %index in entire series
    if t(i) <= mt
        %breaks up bins based on bin length or 15min interval cutoffs
        if t(i) < t0 + bs %&& mod(t(i), bs) %mod breaks up bins on even intervals, non 0 means continue
            ind = ind +1;
            bins{b}.px(ind) = tr.px(i);
            bins{b}.py(ind) = tr.py(i);
            bins{b}.pt(ind) = tr.pt(i);
        else
            ind = 1;
            t0 = t(i);
            b = b + 1;
            bins{b}.px(ind) = tr.px(i);
            bins{b}.py(ind) = tr.py(i);
            bins{b}.pt(ind) = tr.pt(i);
        end
    end
end
%spikes
b = 1; ind = 0;
for i = 1:length(st)
    if st(i) <= mt && st(i) <= max(bins{length(bins)}.pt) %cuts off spikes past max trajectory time in last bin
        while st(i) > max(bins{b}.pt) %skips past bins with no spikes
            b = b + 1;
            ind = 0;
        end
        ind = ind +1;
        bins{b}.sx(ind) = tr.sx(i);
        bins{b}.sy(ind) = tr.sy(i);
        bins{b}.st(ind) = tr.st(i);
    end
end
%remove weak bins less than minspikes, min lengths 
b = 1;
t = {};
for i = 1:length(bins)
    if ~isfield(bins{b},'sx') %% add dummy bins to ends, 297 in middle
        bins{b}.sx(1) = 0; bins{b}.sy(1) = 0; bins{b}.st(1) = bins{b}.pt(1); %adds dummy spike in empty bins
    end
    if isfield(bins{i},'sx') == 1 && length(bins{i}.st) >= minSpikes &&... %if isfield(bins{i},'sx') == 0
            floor((bins{i}.pt(end) - bins{i}.pt(1))/60) >= minBinLength %make bins at least this long
        t{b} = bins{i};
        b = b + 1;
    end
end
bins = t;
for i = 1:length(bins)
    bins{i}.px = bins{i}.px'; bins{i}.py = bins{i}.py'; bins{i}.pt = bins{i}.pt';
    bins{i}.sx = bins{i}.sx'; bins{i}.sy = bins{i}.sy'; bins{i}.st = bins{i}.st';
end
end

function pairs_by_file_index = find_pairs(cells)
pairs_by_file_index = [];
grid_cells = 0;
for i=1:length(cells) %a = load(strcat(params.dir_load,files(i).name));
    a = cells{i};
    %if (strcmp(a.cell_type_obj_new, 'GRID') || strcmp(a.cell_type_obj_new, 'CONJ')) % '_GRID' ?????
    if a.before.gridscore > 0.3
        for j=i+1:length(cells)
            b = cells{j};
            %adjust for missing third <<<<<<
            if (       b.before.gridscore > 0.3 ...
                    && strcmp(a.id, b.id) && strcmp(a.date, b.date)...
                    ... && (strcmp(a.cell_type_subj, 'GRID') || strcmp(a.cell_type_subj, 'CONJ'))  ...
                    && strcmp(a.a, b.a)  ...  %area????
                    && ( a.tet ~= b.tet || a.cell ~= b.cell) ...
                    ... && length(a.B) == 3 &&  length(b.B) == 3 ... %ignore cells without before/after
                    )
                %fprintf('same rat:%s date:%s a:%d,%d b:%d,%d \n', a.rat, ...
                %a.date, a.tetrode, a.cell, b.tetrode, b.cell);
                
                grid_cells = grid_cells + 1;
                pairs_by_file_index(grid_cells,1) = i;
                pairs_by_file_index(grid_cells,2) = j;
            end
        end
    else
        %disp (a.cell_type_subj);
    end
end
fprintf('%d pairs found\n',grid_cells);
end

function pairs_by_file_index = find_pairs_of_cells_mosimol(data)
pairs_by_file_index = [];
grid_cells = 0;
for i=1:length(data) %a = load(strcat(params.dir_load,files(i).name));
    a = data(i); a = a{1}.db;
    if (strcmp(a.cell_type_obj_new, 'GRID') || strcmp(a.cell_type_obj_new, 'CONJ')) % '_GRID' ?????
        for j=i+1:length(data)
            b = data(j); b = b{1}.db;
            if length(a.B) ~= 3
                disp('missing 3rd');
            end
            if ( strcmp(a.rat, b.rat) && strcmp(a.date, b.date)...
                    && (strcmp(a.cell_type_subj, 'GRID') || strcmp(a.cell_type_subj, 'CONJ'))  ...
                    && strcmp(a.area, b.area)  ...  %area????
                    && ( a.tetrode ~= b.tetrode || a.cell ~= b.cell) ...
                    && length(a.B) == 3 &&  length(b.B) == 3 ... %ignore cells without before/after
                    )
                %fprintf('same rat:%s date:%s a:%d,%d b:%d,%d \n', a.rat, ...
                %a.date, a.tetrode, a.cell, b.tetrode, b.cell);
                
                grid_cells = grid_cells + 1;
                pairs_by_file_index(grid_cells,1) = i;
                pairs_by_file_index(grid_cells,2) = j;
            end
        end
    else
        %disp (a.cell_type_subj);
    end
end
fprintf('%d paris found\n',grid_cells);
end




