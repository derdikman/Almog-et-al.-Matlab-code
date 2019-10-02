%{
%}
warning off;
%%%
binsize = 45;
%fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_c_midall_gridscore.mat',binsize); first set
fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45);
load(fn);
[groups ~] = findSimultaneouslyRecordedCells(cells);
gridThreshBef = 0.3; gridThreshMid = 0.25;
[bgids, bcids] = groupsByGridThresh(groups, gridThreshBef, gridThreshMid);
gids = []; cids = [];
[gids, cids] = removeOverlappingCells(groups, bgids, bcids, -1, 0.07,[133]);

%% TIME CORR
%iterate through all groups
%correlate all pairs
%smooth at different sized windows
ijcb = []; ijcm = []; em = [1,2,10,25,50,100,500,2000,5000];%corrBefore = {};  corrMid = {};
tic
for j=1:len(pairs)
    j
    %BEFORE
    tic
    bcorsmooth = []; mcorsmooth = [];
    c1 = cells{pairs(j,1)}; c2 = cells{pairs(j,2)};
    [t1, t2] = createMsSpikeTrain(c1.before.st, 100, c2.before.st);
    for i = 1:len(em)
        t1smooth = movmean(t1,em(i));
        t2smooth = movmean(t2,em(i));
        bcorsmooth(end+1) = corr(t1smooth', t2smooth'); %xcov(train1smooth, train2smooth,0, 'coef');
    end
    %DURING
    [t1, t2] = createMsSpikeTrain(c1.midall.st, 100, c2.midall.st);
    for i = 1:len(em)
        t1smooth = movmean(t1,em(i));
        t2smooth = movmean(t2,em(i));
        mcorsmooth(end+1) = corr(t1smooth', t2smooth');%xcov(train1smooth, train2smooth,0, 'coef');
    end
    %gcorrBefore{k,j,ri} = shuffleTimeCorrelations (g(j).before, g(k).before, n, 200, 0.006);
    %gcorrMid{k,j,ri}     = shuffleTimeCorrelations (g(j).midall, g(k).midall, n, 200, 0.006);
    %corrBefore{g(j).ind,g(k).ind} = bcorsmooth;
    %corrMid{g(j).ind,g(k).ind}     = mcorsmooth;
    ijcb(end+1,:) = [c1.ind,c2.ind, bcorsmooth];
    ijcm(end+1,:) = [c1.ind,c2.ind, mcorsmooth];
    toc
end


%TIME CORR FIG
f = figure(22345);
set(f,'Color','w', 'Position', [200 0 1000 1200]);
for i = 1:len(em)
    a = subplot(3,3,i);
    x = ijcb(:,i+2); y = ijcm(:,i+2);
    plot(x,y,'.','markersize',10); 
    hold on;
    title(sprintf('%dms',em(i))); 
    if i==1;xlabel('PRE'); ylabel('DUR');end
    f1 = fit(x,y,'poly1');
    mm = slim(gca); pf = plot(mm,f1(mm)); pf.LineWidth = lfs-1;
    yfit =  f1.p1 * x + f1.p2;
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    text(a,0.1,0.9,sprintf('r = %.2f res = %.2f',ccof(x, y),rsq),'Units','normalized');
    %legend( sprintf('r = %.2f res = %.2f',ccof(x, y),rsq)  ,'Location','north');
    a.XLim = mm; a.YLim = a.XLim; axis(a,'square');
    hold off
end
suptitle('temporal correlations pre vs dur by smoothing window');

%plot([min(x) max(x)], polyval( polyfit( x, y,1), [min(x) max(x)] )); %'linewidth',lfs,'color','b'
%b = x\y;
%plot(x,b*x);
%text(gca,0.8*max(d.x),polyval(polyfit(d.x, d.y, 1),0.8*max(d.x)), sprintf('%.2f %.2f', round(polyfit(d.x, d.y, 1),2)),'Color','b','FontSize',afs);
%text(sprintf('%.1f %.1f', round( polyfit(x, y, 1), 1 ) ) );
%cc = xcorr2(c2b.rm-mean(c2b.rm(:)),c1b.rm-mean(c1b.rm(:))); %reverse order for perspective
%cc = xcorr2(rm2,rm1); %reverse order for perspective
%cc = normxcorr2(rm2,rm1); %reverse order for perspective
%cenr = round(size(cc,2)/2); cenc = round(size(cc,1)/2);
%ccg = imgaussfilt(cc,sigma); %imagesc(ccg); title(ccg(cenr,cenc));
%tb(end+1) = ccg(cenr,cenc);

%% SPACE CORR
ijspacb = []; ijspacm = []; %e = [1e-9:1:11]%e = 12.^(-1:0.2:1); e = [1e-9 e];
em = [1,2,10,25,50,100,500,2000,5000]; nb = 100;
for i = 1:len(pairs);
    i
    c1 = cells{pairs(i,1)};c2 = cells{pairs(i,2)};
    c1b = cells{pairs(i,1)}.before; c2b = cells{pairs(i,2)}.before;
    t1b=createMsSpikeTrain(c1b.st);t2b=createMsSpikeTrain(c2b.st);
    c1m = cells{pairs(i,1)}.midall; c2m = cells{pairs(i,2)}.midall;
    t1m=createMsSpikeTrain(c1m.st);t2m=createMsSpikeTrain(c2m.st);
    tb = []; tm = [];
    tic
    for j = 1:length(em);
        trs = movmean(t1b,em(j)); rm1 = createSmoothRateMap(c1b,nb,t1b);
        trs = movmean(t2b,em(j)); rm2 = createSmoothRateMap(c2b,nb,trs);
        tb(end+1) =  ccof(rm2,rm1);%cc(cenr,cenc);
        if c1.ind == 34;
            figure(342);
            subplot(3,3,j);
            %imagesc(imgaussfilt(cc,1)); title(ccg(cenr,cenc));colormap('jet');
            imagesc(normxcorr2(rm2,rm1));title(sprintf('%dms',em(j)));colormap('jet');
            hold on;
            %plot(cenr,cenc,'mo')
        end
        %MIDALL
        trs = movmean(t1m,em(j)); rm1 = createSmoothRateMap(c1m,nb,trs);
        trs = movmean(t2m,em(j)); rm2 = createSmoothRateMap(c2m,nb,trs);
        tm(end+1) = ccof(rm2,rm1);%cc(cenr,cenc);
    end
    ijspacb(end+1,:) = [c1.ind,c2.ind, tb];
    ijspacm(end+1,:) = [c1.ind,c2.ind, tm];
    toc
end

SAVE AND LOOK AT NB 100

%SPACE CORR FIG
f=figure(74);set(f,'Color','w','Position',[200 0 1000 1200]);
for i = 1:len(em)
    a=subplot(3,3,i);
    x = ijspacb(:,2+i); y = ijspacm(:,2+i);
    plot(x,y,'.','markersize',10); 
    hold on;
    title(sprintf('%dms',em(i))); 
    if i==1;xlabel('PRE'); ylabel('DUR');end
    f1 = fit(x,y,'poly1');
    mm = slim(gca); pf = plot(mm,f1(mm)); pf.LineWidth = lfs-1;
    yfit =  f1.p1 * x + f1.p2;
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    text(a,0.1,0.9,sprintf('r = %.2f res = %.2f',ccof(x, y),rsq),'Units','normalized');
    a.XLim = mm; a.YLim = a.XLim; axis(a,'square');
    hold off
end
suptitle('spatial correrlations pre vs dur by smoothing window');



tic
%spike removal
for ri = 1:length(bgids)
    %ri
    %gis = gids(ri); cis = cids{gids(ri)};%gis}
    g = groups{bgids(ri)};%gis};
    g = g(bcids{bgids(ri)});%cis);
    for j = 1:length(g)-1
        for k = j+1:length(g)
            %BEFORE
            fprintf('b %3dx%3d: ',g(j).ind,g(k).ind);
            c1 = g(j).before; c2 = g(k).before;
            spkt1=floor(c1.st*1000)+1; spkt2=floor(c2.st*1000)+1; min_time=min(min(spkt1),min(spkt2)); %TIME IN MS
            offset = 100; spkt1 = spkt1-min_time+offset; spkt2 = spkt2-min_time+offset; max_time=max(max(spkt1),max(spkt2))+offset;
            [s1i, s2i] = removeOverlappingSpikes(spkt1,spkt2, 1);
            c1 = g(j).midall; c2 = g(k).midall;
            spkt1=floor(c1.st*1000)+1; spkt2=floor(c2.st*1000)+1; min_time=min(min(spkt1),min(spkt2)); %TIME IN MS
            offset = 100; spkt1 = spkt1-min_time+offset; spkt2 = spkt2-min_time+offset; max_time=max(max(spkt1),max(spkt2))+offset;
            fprintf('m %3dx%3d: ',g(j).ind,g(k).ind);
            [s1i, s2i] = removeOverlappingSpikes(spkt1,spkt2, 1);
            %disp('_')
        end
    end
end
toc;
%plot(size(cc,2)/2,size(cc,1)/2,'md','MarkerFaceColor','w','MarkerSize',7);
%imagesc(az, imgaussfilt(c1.rm,1)); %RM


% ijcb = []; ijcm = [];
%     for ri = 1:283
%         for ci = 1:287
%             if ~isempty(corrBefore{ri,ci})
%                 t = corrBefore{ri,ci};
%                 ijcb(end+1,:) = [ri,ci,t];
%                 t = corrMid{ri,ci};
%                 ijcm(end+1,:) = [ri,ci,t];
%                 %ijcbmpbm(end+1,:) = [ri ci cb cm pb pm ];
%             end
%         end
%     end



%{

%TIME CORR FIG
f = figure(22345);
set(f,'Color','w', 'Position', [200 0 1200 1200]);
for i = 1:len(em)
    subplot(3,3,i);
    x = ijcb(:,i+2); y = ijcm(:,i+2); a = -99; b = 2; c = length(x);
    %y(x>b*std(x)) = a; x(x>b*std(x)) = a; x(y>b*std(y)) = a; y(y>b*std(y)) = a; x = x(x~=a); y = y(y~=a); %go by std of x
    %fprintf('%dms corrs removed: %d\n',2^(i-1),c - length(x));
    plot(x,y,'.');    %'mo','markersize',16,'linewidth',lfs
    hold on; %axis('tight');s
    p = polyfit(x,y,1);
    %yfit = polyval(p,x);
    yfit =  p(1) * x + p(2);
    plot(x,yfit);
    % %FROM F3
    f1 = fit(ctsbma(:,4),            ctsbma(:,6),'poly1');
    mm = slim(axD); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
    %
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    %legend((sprintf('r = %.1f %.1f', polyfit( x, y,1) )  ),'Location','best');
    %legend( sprintf('a = %.2f res = %.2f',p(1),rsq)  ,'Location','north');
    lim = [min([x;y]); max([x;y])]; xlim(lim); ylim(lim);
    %textbest(gca, sprintf('a = %.2f\nres = %.2f',p(1),rsq),{'Location','bestoutside'},{});
    title(sprintf('%dms r=%.2f res=%.2f',em(i),ccof(x,y),rsq));
    hold off;axis('square');%axis('equal');
end
suptitle('temporal correlations pre vs dur by smoothing window');
%}