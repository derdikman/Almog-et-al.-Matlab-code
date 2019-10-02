%UNCOMMENT TO LOAD CELLS
%load(sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45));

%% preamble
dbstop if error
warning off
loads();
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
    ijcb(end+1,:) = [c1.ind,c2.ind, bcorsmooth];
    ijcm(end+1,:) = [c1.ind,c2.ind, mcorsmooth];
    toc
end


%TIME CORR FIG
f = figure(9331);
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
    plot(x,yfit)
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    %legend((sprintf('r = %.1f %.1f', polyfit( x, y,1) )  ),'Location','best');
    %legend( sprintf('a = %.2f res = %.2f',p(1),rsq)  ,'Location','north');
    lim = [min([x;y]); max([x;y])]; xlim(lim); ylim(lim);
    %textbest(gca, sprintf('a = %.2f\nres = %.2f',p(1),rsq),{'Location','bestoutside'},{});
    title(sprintf('%dms a=%.2f res=%.2f',em(i),p(1),rsq));
    hold off;axis('square');%axis('equal');
end
suptitle('time corr bef vs dur by smoothing');

%plot([min(x) max(x)], polyval( polyfit( x, y,1), [min(x) max(x)] )); %'linewidth',lfs,'color','b'
%b = x\y;
%plot(x,b*x);
%text(gca,0.8*max(d.x),polyval(polyfit(d.x, d.y, 1),0.8*max(d.x)), sprintf('%.2f %.2f', round(polyfit(d.x, d.y, 1),2)),'Color','b','FontSize',afs);
%text(sprintf('%.1f %.1f', round( polyfit(x, y, 1), 1 ) ) );


%SPACE CORR
%do spatial smoothing, sample middle do before and after
ijspacb = []; ijspacm = []; %e = [1e-9:1:11]%e = 12.^(-1:0.2:1); e = [1e-9 e];
em = [1,2,10,25,50,100,500,2000,10000];
for i = 1:len(pairs);
    i
    c1 = cells{pairs(i,1)};c2 = cells{pairs(i,2)};
    %[t1m, t2m] = createMsSpikeTrain(c1m.st, 100, c2m.st); %NO OFFS
    c1b = cells{pairs(i,1)}.before; c2b = cells{pairs(i,2)}.before;
    t1b=createMsSpikeTrain(c1b.st);t2b=createMsSpikeTrain(c2b.st);
    c1m = cells{pairs(i,1)}.midall; c2m = cells{pairs(i,2)}.midall;
    t1m=createMsSpikeTrain(c1m.st);t2m=createMsSpikeTrain(c2m.st);
    tb = []; tm = [];
    tic
    for j = 1:length(em);
        %sigma = em(i);
        %c1 = pairs(i,1).before; c2 = pairs(i,2).before;
        trs = movmean(t1b,em(j)); rm1 = createSmoothRateMap(c1b,50,t1b);
        trs = movmean(t2b,em(j)); rm2 = createSmoothRateMap(c2b,50,trs);
        %cc = xcorr2(c2b.rm-mean(c2b.rm(:)),c1b.rm-mean(c1b.rm(:))); %reverse order for perspective
        %cc = xcorr2(rm2,rm1); %reverse order for perspective
        cc = normxcorr2(rm2,rm1); %reverse order for perspective
        cenr = round(size(cc,2)/2); cenc = round(size(cc,1)/2);
        %ccg = imgaussfilt(cc,sigma); %imagesc(ccg); title(ccg(cenr,cenc));
        %tb(end+1) = ccg(cenr,cenc);
        tb(end+1) = cc(cenr,cenc);
        if c1.ind == 34;
            figure(9933);
            subplot(3,3,j);
            %imagesc(imgaussfilt(cc,1)); title(ccg(cenr,cenc));colormap('jet');
            imagesc(cc);title(sprintf('%dms',em(j)));colormap('jet');
            hold on;
            plot(cenr,cenc,'mo')
            
        end
        %imagesc(cc); %xlabel(az, sprintf('c%d x c%d',c1.ind,c2.ind),'fontweight','bold');colormap(az,'jet');
        %imagesc(imgaussfilt(cc,1)); title(ccg(cenr,cenc))%colormap('jet');
        %MIDALL
        %c1 = pairs(i,1).midall; c2 = pairs(i,2).midall;
        %cc = xcorr2(c2.rm-mean(c2.rm(:)),c1.rm-mean(c1.rm(:))); %reverse order for perspective
        trs = movmean(t1m,em(j)); rm1 = createSmoothRateMap(c1m,50,trs);
        trs = movmean(t2m,em(j)); rm2 = createSmoothRateMap(c2m,50,trs);
        %cc = xcorr2(c2.rm-mean(c2.rm(:)),c1.rm-mean(c1.rm(:))); %reverse order for perspective
        cc = normxcorr2(rm2,rm1); %reverse order for perspective
        cenr = round(size(cc,2)/2); cenc = round(size(cc,1)/2);
        %ccg = imgaussfilt(cc,sigma);
        %tm(end+1) = ccg(cenr,cenc);
        tm(end+1) = cc(cenr,cenc);
    end
    ijspacb(end+1,:) = [c1.ind,c2.ind, tb];
    ijspacm(end+1,:) = [c1.ind,c2.ind, tm];
    toc
end


%SPACE CORR FIG
figure(9932);
for i = 1:len(em)
    subplot(3,3,i);
    sigma = em(i);
    %subplot(3,6,i);
    x = ijspacb(:,2+i); y = ijspacm(:,2+i); a = 0; b = 1.6;
    %y(x>b*std(x)) = a; x(x>b*std(x)) = a; x(y>b*std(y)) = a; y(y>b*std(y)) = a; x = x(x~=0); y = y(y~=0); %go by std of x
    plot(x,y,'.'); axis('tight');a = xlim; b = ylim; xlim([min(a(1),b(1)), max(a(2),b(2))]); ylim(xlim)   %'mo','markersize',16,'linewidth',lfs
    hold on;
    title(sprintf('%dms',em(i))); xlabel('bef'); ylabel('during');
    p = polyfit(x,y,1);
    yfit = polyval(p,x);yfit =  p(1) * x + p(2);
    plot(x,yfit)
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    %legend((sprintf('r = %.1f %.1f', polyfit( x, y,1) )  ),'Location','best');
    legend( sprintf('a = %.2f res = %.2f',p(1),rsq)  ,'Location','best');
    hold off;%axis('equal');
    axis('square');
    
end
suptitle('spatial corr at 0,0 bef vs dur by spike train smoothing window');


