%% grid thresh
t = [];
for j=1:len(cells)
    c = cells{j}; p=[]; p.movmean = 25; nbins = 50; %spatial    
    gd(j,:) = [j c.before.gridscore c.midall.gridscore];
    t{j} = [c.id c.date];
end
t1 = 0.5; t2 = 0.0; %0.3 0.25
tt = gd(:,2)>=t1 & gd(:,3)<=t2;
celst=gd(tt,1); 
cc = t(tt);
pairst = [];
for i = 1:len(celst)
    for j = i+1:len(celst)
        if strcmp(cc{i}, cc{j})
            pairst(end+1,:) = [celst(i) celst(j)];
        end
    end
end
pairst = removeOverlappingPairs(cells, pairst, 1, 0.7);
[~,t,~] = intersect(pairs,pairst,'rows'); %which pairs to use ..

ccofts = [ccof(ctsbma(t,1),ctsbma(t,2)), ccof(ctb(t), ctm(t)),... dfdfsdf
ccof(ctsbma(t,4),ctsbma(t,5)), ccof(csb(t),csm(t))]

ccofts = [ccof(ctsbma(:,1),ctsbma(:,2)), ccof(ctb, ctm),... dfdfsdf
ccof(ctsbma(:,4),ctsbma(:,5)), ccof(csb,csm)]

figure(665); clf; 
subplot(221);d1 = ctsbma(t,1); d2 = ctsbma(t,2); plot(d1,d2,'.'); axis('tight'); 
title(sprintf('grid thresh pre >= %.2f dur <=%.2f\ntime',t1,t2));xlabel('pre'); ylabel('dur'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-'); axis('tight'); axis square;
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(223);d1 = ctsbma(t,4); d2 = ctsbma(t,5); plot(d1,d2,'.'); axis('tight'); 
title('space');xlabel('pre'); ylabel('dur'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');  axis('tight'); axis square;
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(222);d1 = ctsbma(:,1); d2 = ctsbma(:,2); plot(d1,d2,'.'); axis('tight'); 
title(sprintf('grid thresh pre >= %.2f dur <=%.2f\ntime',0.3,0.25));xlabel('pre'); ylabel('dur'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-'); axis('tight'); axis square;
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(224);d1 = ctsbma(:,4); d2 = ctsbma(:,5); plot(d1,d2,'.'); axis('tight'); 
title('space');xlabel('pre'); ylabel('dur'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-'); axis('tight'); axis square;
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

suptitle('Changing Grid Thresholds, col1 new, col2 original');



%% grid
tic; gd=[];
for j=1:len(cels)
    c = cells{cels(j)}; p=[]; p.movmean = 25; nbins = 50; %spatial
    strt = 45 *60 -inf; midt = (45/1) *60; endt=45 *60 +inf;
    ca = windowsesh(c.midall,-1,-1,strt,midt);
    cb = windowsesh(c.midall,-1,-1,midt,endt);
    %rm = createRateMap(s.px, s.py, s.pt, s.sx, s.sy, s.st, true, n);
    %s.ac = xcorr2(s.rm);%s.gridscore = gridscore2(s.ac, 2);
    gd(j,:) = [gridscore2(xcorr2(createRateMap(ca,nbins)),2),gridscore2(xcorr2(createRateMap(cb,nbins)),2)];    
end
toc

%% WINDOWS
tic; wd=[]; wb=[]; %was 100 (anylysis done in 50?)
for j=1:len(pairs)
    j
    c1 = cells{pairs(j,1)}; c2 = cells{pairs(j,2)}; p=[]; p.movmean = 25; nbins = 50; %spatial
%     wb(j) = c1.midall.pt(end)-c1.midall.pt(1);
%    t = c1.midall.pt; tm = t(1)+(t(end)-t(1))/2;
%     [train1,train2] = createMsSpikeTrain(c1.midall.st,t(end),c2.midall.st, t(1),tm);
%     [a1, ~, t1sm, t2sm] = timeCorrelationSmoothed( train1,train2,p);
    
    %im = ceil(len(c1.midall.pt)/2); find index of half way
    strt = 45 *60 -inf; midt = (45/1) *60; endt=45 *60 +inf;
    %is = find(c1.midall.pt>=strt*60,1);
    %im = find(c1.midall.pt>=midt*60,1);
    %ie = len(c1.midall.pt);
    %ie = find(c1.midall.pt<=endt*60,1,'last');
    
    
    %c1a = windowsesh(c1.midall,is,im); c2a = windowsesh(c2.midall,is,im); %was 1
    %c1b = windowsesh(c1.midall,im,ie); c2b = windowsesh(c2.midall,im,ie);
    c1a = windowsesh(c1.midall,0,0,strt,midt); c2a = windowsesh(c2.midall,0,0,strt,midt); %was 1
    c1b = windowsesh(c1.midall,0,0,midt,endt); c2b = windowsesh(c2.midall,0,0,midt,endt);
    
    if ~(isempty(c1a) || isempty(c1b) || isempty(c2a) || isempty(c2b) )
        %1st half
        [t1,t2] = createMsSpikeTrain(c1a.st,-1,c2a.st);
        [ta, ~, t1sm, t2sm] = timeCorrelationSmoothed(t1,t2,p);
        sa=ccof(createSmoothRateMap(c2a,nbins,t2sm),createSmoothRateMap(c1a,nbins,t1sm));
        %2nd half
        [t1,t2] = createMsSpikeTrain(c1b.st,-1,c2b.st);
        [tb, ~, t1sm, t2sm] = timeCorrelationSmoothed(t1,t2,p);
        sb=ccof(createSmoothRateMap(c2b,nbins,t2sm),createSmoothRateMap(c1b,nbins,t1sm));
    else
        ta=0;tb=0;sa=0;sb=0;
    end
    sa(sa==1)=0;sb(sb==1)=0;
    wd(j,:) = [ta tb sa sb];  
end
wd(isnan(wd))=0; toc

figure(66353); clf; xl = 'corr pre'; yl = 'corr dur';
%TIME
%1st half
subplot(221);d1 = ctsbma(:,1); d2 = wd(:,1); plot(d1,d2,'.'); axis('tight'); axis square;
title(sprintf('[%.0f : %.0f min]',strt/60,midt/60));
xlabel(xl); ylabel(yl); hold on; 
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
%2nd half
subplot(222);d1 = ctsbma(:,1); d2 = wd(:,2); plot(d1,d2,'.'); axis('tight'); axis square;
title(sprintf('[%.0f : %.0f min]',midt/60,endt/60));
xlabel(xl); ylabel(yl); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
%SPACE
subplot(223);d1 = ctsbma(:,4); d2 = wd(:,3); plot(d1,d2,'.'); axis('tight'); axis square;
title(' ');xlabel(xl); ylabel(yl); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
%2nd half
subplot(224);d1 = ctsbma(:,4); d2 = wd(:,4); plot(d1,d2,'.'); axis('tight'); axis square;
title(' ');xlabel(xl); ylabel(yl); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

suptitle('correlations by muscimol time window (row 1 temporal, row 2 spatial)');



%% MEAN FIRING RATE
fb=[];fd=[];fa=[];
for j=1:len(pairs)
    c1 = cells{pairs(j,1)}; c2 = cells{pairs(j,2)};
    s1 = c1.before.st; s1= len(s1)/(s1(end)-s1(1)); 
    s2 = c2.before.st; s2= len(s2)/(s2(end)-s2(1)); 
    fb(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
    s1 = c1.midall.st; s1= len(s1)/(s1(end)-s1(1)); 
    s2 = c2.midall.st; s2= len(s2)/(s2(end)-s2(1));     
    fd(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
    s1 = c1.after.st; s1= len(s1)/(s1(end)-s1(1)); 
    s2 = c2.after.st; s2= len(s2)/(s2(end)-s2(1));
    fa(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
end
fa(fa(:,1)==inf,:)=[];

rsp=[];
for i=1:length(pairs)
    rsp(i,:) = [cells{pairs(i,1)}.midall.rayleigh_score,...
                cells{pairs(i,2)}.midall.rayleigh_score];
end
clhd = rsp(:,1)>0.4 & rsp(:,2)>0.4;
%FIRING RATE FIG

figure(662); clf; xl = 'mean rate (Hz)'; yl = 'corr';
%BEFORE
% g score
subplot(231);d1 = fb(:,3); d2 = mb(:,3); plot(d1,d2,'.'); axis('tight');
title('PRE:: grid score');xlabel(xl); ylabel('score'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'k-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
[];plot(d1(clhd), d2(clhd),'r.'); f1 = fit(d1(clhd), d2(clhd),'poly1'); plot(d1,f1(d1),'r-');
text(0.1,0.8,sprintf('hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(clhd),d2(clhd)),3)),...
    'Units','normalized','color','r');
f1 = fit(d1(~clhd), d2(~clhd),'poly1'); plot(d1,f1(d1),'b-');
text(0.1,0.7,sprintf('~hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(~clhd),d2(~clhd)),3)),...
    'Units','normalized','color','b');
%time
subplot(232);d1 = fb(:,3); d2 = ctsbma(:,1); plot(d1,d2,'.'); axis('tight'); 
title('temporal');xlabel(xl); ylabel(yl); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'k-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
[];plot(d1(clhd), d2(clhd),'r.'); f1 = fit(d1(clhd), d2(clhd),'poly1'); plot(d1,f1(d1),'r-');
text(0.1,0.8,sprintf('hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(clhd),d2(clhd)),3)),...
    'Units','normalized','color','r');
f1 = fit(d1(~clhd), d2(~clhd),'poly1'); plot(d1,f1(d1),'b-');
text(0.1,0.7,sprintf('~hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(~clhd),d2(~clhd)),3)),...
    'Units','normalized','color','b');
%space
subplot(233);d1 = fb(:,3); d2 = ctsbma(:,4); plot(d1,d2,'.'); axis('tight'); 
title('spatial');xlabel(xl); ylabel(yl); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'k-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
[];plot(d1(clhd), d2(clhd),'r.'); f1 = fit(d1(clhd), d2(clhd),'poly1'); plot(d1,f1(d1),'r-');
text(0.1,0.8,sprintf('hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(clhd),d2(clhd)),3)),...
    'Units','normalized','color','r');
f1 = fit(d1(~clhd), d2(~clhd),'poly1'); plot(d1,f1(d1),'b-');
text(0.1,0.7,sprintf('~hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(~clhd),d2(~clhd)),3)),...
    'Units','normalized','color','b');
%DURING
%score
subplot(234);d1 = fd(:,3); d2 = md(:,3); plot(d1,d2,'.'); axis('tight'); 
title('DUR:: grid score');xlabel(xl); ylabel('score'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'k-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
[];plot(d1(clhd), d2(clhd),'r.'); f1 = fit(d1(clhd), d2(clhd),'poly1'); plot(d1,f1(d1),'r-');
text(0.1,0.8,sprintf('hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(clhd),d2(clhd)),3)),...
    'Units','normalized','color','r');
f1 = fit(d1(~clhd), d2(~clhd),'poly1'); plot(d1,f1(d1),'b-');
text(0.1,0.7,sprintf('~hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(~clhd),d2(~clhd)),3)),...
    'Units','normalized','color','b');
%time
subplot(235);d1 = fd(:,3); d2 = ctsbma(:,2); plot(d1,d2,'.'); axis('tight'); 
title('temporal');xlabel(xl); ylabel(yl); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'k-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
[];plot(d1(clhd), d2(clhd),'r.'); f1 = fit(d1(clhd), d2(clhd),'poly1'); plot(d1,f1(d1),'r-');
text(0.1,0.8,sprintf('hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(clhd),d2(clhd)),3)),...
    'Units','normalized','color','r');
f1 = fit(d1(~clhd), d2(~clhd),'poly1'); plot(d1,f1(d1),'b-');
text(0.1,0.7,sprintf('~hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(~clhd),d2(~clhd)),3)),...
    'Units','normalized','color','b');
%space
subplot(236);d1 = fd(:,3); d2 = ctsbma(:,5); plot(d1,d2,'.'); axis('tight'); 
title('spatial');xlabel('mean rate'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'k-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
[];plot(d1(clhd), d2(clhd),'r.'); f1 = fit(d1(clhd), d2(clhd),'poly1'); plot(d1,f1(d1),'r-');
text(0.1,0.8,sprintf('hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(clhd),d2(clhd)),3)),...
    'Units','normalized','color','r');
f1 = fit(d1(~clhd), d2(~clhd),'poly1'); plot(d1,f1(d1),'b-');
text(0.1,0.7,sprintf('~hd: a=%.3f r=%.3f',round(f1.p1,3),round(ccof(d1(~clhd),d2(~clhd)),3)),...
    'Units','normalized','color','b');

suptitle('Mean Firing Rate (row 1 before, row 2 during)');

a = [];
for j=1:len(cels)
    c = cells{j}.before;
    a(j,:) = [1/(min(diff(c.st))),c.max_r]; %max rate
end

for j=1:len(cels)
c = cells{cels(j)}.before; b=1./diff(c.st);
[sum(b<1) sum(b>1) sum(b>100) sum(b>500) mean(b(b<500 & b>1))]
end;



%% SCORE VS CORR 
mb=[];md=[];
for j=1:len(pairs)
    c1 = cells{pairs(j,1)}; c2 = cells{pairs(j,2)};%c1.ind, c2.ind
    s1 = c1.before.gridscore; s2 = c2.before.gridscore; s1(s1==-2)=0; s2(s2==-2)=0;
    mb(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
    s1 = c1.midall.gridscore; s2 = c2.midall.gridscore; s1(s1==-2)=0; s2(s2==-2)=0;
    md(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
end


figure(661); clf; tos=4;
subplot(231);d1 = mb(:,1); d2 = ctsbma(:,tos); plot(d1,d2,'.'); axis('tight'); 
title('min gs');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-'); 
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(232);d1 = mb(:,2); d2 = ctsbma(:,tos); plot(d1,d2,'.'); axis('tight'); 
title('max gs');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1'); plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(233);d1 = mb(:,3); d2 = ctsbma(:,tos); plot(d1,d2,'.'); axis('tight'); 
title('mean gs');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(234);d1 = md(:,1); d2 = ctsbma(:,tos+1); plot(d1,d2,'.'); axis('tight'); 
title('min gs');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(235);d1 = md(:,2); d2 = ctsbma(:,tos+1); plot(d1,d2,'.'); axis('tight'); 
title('max gs');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(236);d1 = md(:,3); d2 = ctsbma(:,tos+1); plot(d1,d2,'.'); axis('tight'); 
title('mean gs');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

if tos==1
suptitle('Grid Score vs Time Corr (row 1 before, row 2 during)');
else
    suptitle('Grid Score vs Space Corr (row 1 before, row 2 during)');
end

%% SUPP SCORE VS CORR 
mb=[];md=[];
for j=1:len(pairs)
    c1 = cells{pairs(j,1)}; c2 = cells{pairs(j,2)};%c1.ind, c2.ind
    s1 = c1.before.gridscore; s2 = c2.before.gridscore; s1(s1==-2)=0; s2(s2==-2)=0;
    mb(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
    s1 = c1.midall.gridscore; s2 = c2.midall.gridscore; s1(s1==-2)=0; s2(s2==-2)=0;
    md(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
end


figure(6611); clf; tos=1;

subplot(241);d1 = mb(:,3); d2 = ctsbma(:,tos); plot(d1,d2,'.'); axis('tight'); 
title('time corr pre');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(242);d1 = md(:,3); d2 = ctsbma(:,tos+1); plot(d1,d2,'.'); axis('tight'); 
title('time corr dur');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

tos = 4;

subplot(243);d1 = mb(:,3); d2 = ctsbma(:,tos); plot(d1,d2,'.'); axis('tight'); 
title('space corr pre');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(244);d1 = md(:,3); d2 = ctsbma(:,tos+1); plot(d1,d2,'.'); axis('tight'); 
title('space corr dur');xlabel('score'); ylabel('corr'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

tos=1;

subplot(245);d1 = mb(:,3); d2 = ctsbma(:,tos); d2 = abs(d2); plot(d1,d2,'.'); axis('tight'); 
title('time corr pre');xlabel('score'); ylabel('abs(corr)'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(246);d1 = md(:,3); d2 = ctsbma(:,tos+1); d2 = abs(d2); plot(d1,d2,'.'); axis('tight'); 
title('time corr dur');xlabel('score'); ylabel('abs(corr)'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

tos = 4;

subplot(247);d1 = mb(:,3); d2 = ctsbma(:,tos); d2 = abs(d2); plot(d1,d2,'.'); axis('tight'); 
title('space core pre');xlabel('score'); ylabel('abs(corr)'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

subplot(248);d1 = md(:,3); d2 = ctsbma(:,tos+1); d2 = abs(d2); plot(d1,d2,'.'); axis('tight'); 
title('space core dur');xlabel('score'); ylabel('abs(corr)'); hold on;
f1 = fit(d1,d2,'poly1');plot(d1,f1(d1),'-');
text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');

suptitle('Mean Grid Score vs Correlation');
