% uncomment to load cells
%load('Z:\\data\\noam\\Data\\muscimol\\noam\\cells_45min_d_patchtraj_rayleigh'); %lab
%load('C:\\Noam\\Data\\muscimol\\noam\\cells_45min_d_patchtraj_rayleigh'); %personal
%%preamble 
%{ 
 load('C:\Noam\Data\muscimol\cells15nan');
 bthresh=0.5; mthresh=0.2; %MAKE SURE TO MATCH PAIRS MAKE
 [pairs,cels,group,chd,phd,chdn,phdn]=pairsMake(cellsn, bthresh, mthresh)
%}
function f1(cellsn, cels, bthresh, mthresh)
dbstop if error  

fig = figure(991);  fs = 12;
set(fig,'Color','w', 'Position', [600 0 1400 800]);
gtop = uix.GridFlex('Parent',fig,'Spacing',5, 'BackgroundColor','w');
gl = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gl,'BackgroundColor','w'),'visible','off');
text(0,0.5,'A','fontweight','bold','fontsize',fs);
gA = uix.GridFlex('Parent',gl,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');


set(gl,'Widths', [-1], 'Heights', [ 25, -1]);
gr = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');

%1A - trajectory / auto for group
ii = [100, 103, 104, 105, 111];
axg = []; ug = []; sb = []; sm = sb;
sb.x='position(cm)'; sb.y = sb.x;
sm.x=''; sm.y = sm.x;
t=0;
for i = 1:len(ii)
     c = cellsn(ii(i));
     t=max([t max(c.before.px) max(c.midall.px) max(c.after.px)]);
     t=max([t max(c.before.py) max(c.midall.py) max(c.after.py)]);
end
sb.lim = t;sm.lim = t;
for i = 1:len(ii)
    c = cellsn(ii(i));
    sb.t = ' '; 
    axg(end+1) = axes('Parent',uicontainer('Parent',gA)); if(i==3); sb.t = 'pre'; end 
    axg(end)=plotTR(axg(end),c.before,sb);
    sb.x=''; sb.y = sb.x;
    axg(end+1) = axes('Parent',uicontainer('Parent',gA));  if(i==3); sb.t = 'dur'; end 
    axg(end)=plotTR(axg(end),c.midall,sb);
    axg(end+1) = axes('Parent',uicontainer('Parent',gA));  if(i==3); sb.t = 'post'; end 
    axg(end)=plotTR(axg(end),c.after,sb);
    g = [c.before.gridscore c.midall.gridscore c.after.gridscore]; g(g==-2)=0;
    axg(end+1) = axes('Parent',uicontainer('Parent',gA));  sm.t = sprintf('%s%.1f','grid=',g(1));
    axg(end)=plotAC(axg(end),c.before,sm);
    axg(end+1) = axes('Parent',uicontainer('Parent',gA));  sm.t = sprintf('%s%.1f','grid=',g(2));
    axg(end)=plotAC(axg(end),c.midall,sm);
    axg(end+1) = axes('Parent',uicontainer('Parent',gA));  sm.t = sprintf('%s%.1f','grid=',g(3));
    axg(end)=plotAC(axg(end),c.after,sm);
end
set(gA,'Widths', [-1 -1 -1 -1 -1], 'Heights', [-1 -1 -1 -1 -1 -1])

%1B - mean firing pre dur
axes('Parent',uicontainer('Parent',gr,'BackgroundColor','w'),'visible','off');
text(0,0.5,'B','fontweight','bold','fontsize',fs);
uB = uicontainer('Parent',gr); axB = axes(uB); mrbda = [];
for i = 1:len(cels)
    c = cellsn(cels(i));
    mrbda(i,1) = len(c.before.st)/(c.before.st(end)-c.before.st(1));
    mrbda(i,2) = len(c.midall.st)/(c.midall.st(end)-c.midall.st(1));
    mrbda(i,3) = len(c.after.st)/(c.after.st(end)-c.after.st(1));
end
d1 = mrbda(:,1); d2=mrbda(:,2);
plot(d1,d2,'ko'); hold on;
% plot(mrbda(:,3),mrbda(:,2),'bo');
axB.XLim = slim(axB); axB.YLim = axB.XLim; axis(axB,'square');
f = fit(d1,d2,'poly1');plot(axB.XLim,f(axB.XLim),'-'); [r p]=ccof(d1,d2);
text(0.1,0.9,sprintf('a=%.2f r=%.2f %s',rnd(f.p1,2), rnd(r,2), pstr(p,2)),'Units','normalized');
%text(0.1,0.9,sprintf('a=%.2f r=%.2f p<%.3f',rnd(f.p1,2), rnd(r,2), rnd(p,2)),'Units','normalized');
%text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
hold off
xlabel('pre (Hz)'); ylabel('dur (Hz)');
title('mean firing rate pre vs during');
%legend(axB, {'pre','post'});


%1D threshold all
axes('Parent',uicontainer('Parent',gr,'BackgroundColor','w'),'visible','off');
text(0,0.5,'D','fontweight','bold','fontsize',fs);
uD = uicontainer('Parent',gr); axD = axes(uD);
sigbma = [];
for i = 1:len(cels)
    c = cellsn(cels(i));
    sigbma(i,:) = [c.before.gridscore c.midall.gridscore c.after.gridscore]; %?????  
end
agbm = [];
for i = 1:len(cellsn)
    c = cellsn(i);
    agbm(i,:) = [c.before.gridscore c.midall.gridscore];     
end
plot(agbm(:,1),agbm(:,2),'b.'); hold on
plot(sigbma(:,1),sigbma(:,2),'ro');
%plot(gbm((gbm(:,2) == 0),1),gbm((gbm(:,2) == 0),2),'rx');
l= slimd(axD); axD.XLim = l; axD.YLim = l;
plot([bthresh bthresh],[l(1) mthresh],'g');plot([bthresh l(2)],[mthresh mthresh],'g');
hold off; axis(axD,'square');
xlabel('gridscore pre'); ylabel('gridscore dur');
title('gridscore all vs cohort'); 
%legend(axC, {'all','significant','thresholds'});

%1C - mean firing pre post
axes('Parent',uicontainer('Parent',gr,'BackgroundColor','w'),'visible','off');
text(0,0.5,'C','fontweight','bold','fontsize',fs);
uC = uicontainer('Parent',gr); axC = axes(uC);
d1=mrbda(:,1); d2=mrbda(:,3); t=mrbda(:,3)~=inf; d1=d1(t);d2=d2(t);%d1(d1==inf)=0;d2(d2==inf)=0;
plot(d1,d2,'ko'); hold on;
axC.XLim = slim(axC); axC.YLim = axC.XLim; axis(axC,'square');
f = fit(d1,d2,'poly1');plot(axC.XLim,f(axC.XLim),'-');[r p]=ccof(d1,d2);
text(0.1,0.9,sprintf('a=%.2f r=%.2f %s',rnd(f.p1,2), rnd(r,2), pstr(p,2)),'Units','normalized');
%text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
hold off
xlabel('pre (Hz)'); ylabel('post (Hz)');
title('mean firing rate pre vs post');

%1E Significance cels only
axes('Parent',uicontainer('Parent',gr,'BackgroundColor','w'),'visible','off');
text(0,0.5,'E','fontweight','bold','fontsize',fs);
uE = uicontainer('Parent',gr,'BackgroundColor','w');axE = axes(uE); 
d1=sigbma(:,1);d2=sigbma(:,3); t=mrbda(:,3)~=inf; d1=d1(t);d2=d2(t);%d1(d1==inf)=0;d2(d2==inf)=0;
plot(d1,d2,'ko'); hold on;
axE.XLim = slim(axE); axE.YLim = axE.XLim; axis(axE,'square');
f = fit(d1,d2,'poly1');plot(axE.XLim,f(axE.XLim),'-');[r p]=ccof(d1,d2);
text(0.1,0.9,sprintf('a=%.2f r=%.2f %s',rnd(f.p1,2), rnd(r,2), pstr(p,2)),'Units','normalized');
%text(0.1,0.9,sprintf('a=%.2f r=%.2f p=%.2f',rnd(f.p1,2), rnd(r,2), rnd(p,2)),'Units','normalized');
%text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
hold off
xlabel('gridscore pre'); ylabel('gridscore post');
title('gridscore pre vs post');
set( gr, 'Heights', [25 -1 25 -1] );
a = findobj(fig,'type','UIContainer');
for i = 1:len(a)
    a(i).BackgroundColor = 'w';
end

end
%print(f, 'f1.png', '-dpng', '-r300' );   %save file as PNG w/ 300dpi

