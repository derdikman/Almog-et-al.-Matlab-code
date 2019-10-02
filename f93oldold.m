%% preamble
dbstop if error  
%load('.\\data\\figdata');
f = figure(9931); tic
set(f,'Color','w', 'Position', [200 0 1000 1000]);
fs = 14; afs = 10; mfs = 10; tfs = 12; lfs = 3;  p = {};
gtop = uix.GridFlex('Parent', f,'Spacing',5, 'BackgroundColor','w','DividerMarkings','on');

%{


%Pre: 18, 49, 0, 32 space
%Dur:  7, 31, 2, 59

'sig:sig   sig:non   non:sig   non:non'
[sum(pptsbma(:,a)&pptsbma(:,b)),...
sum(pptsbma(:,a)&~pptsbma(:,b)),...
sum(~pptsbma(:,a)&pptsbma(:,b)),...
sum(~pptsbma(:,a)&~pptsbma(:,b)) ]
if(a==1) disp('pre'); else disp('dur'); end

%Pre=  bar([40 27 17 15; 1 37 0 61]');
%Dur=  bar([1, 37, 0, 61])
% title('significance [temp:xgs]')
%xlabel(['sig:sig   ', 'sig:non   ', 'non:sig   ', 'non:non']);
%legend({'PRE','DUR'})
figure(9932);
a=1; b=4; yl = 'space'; t = 'PRE'
a=2; b=5; yl = 'space'; t = 'DUR'; figure(99322);
%  plot(ctsbma(:,a),            ctsbma(:,b),'.','markersize',mfs,'linewidth',lfs); hold on;
plot(ctsbma(pptsbma(:,a)&pptsbma(:,b),a), ctsbma(pptsbma(:,a)&pptsbma(:,b),b),'o','markersize',mfs-3.2,'linewidth',lfs-1.5); hold on;
plot(ctsbma(pptsbma(:,a)&~pptsbma(:,b),a), ctsbma(pptsbma(:,a)&~pptsbma(:,b),b),'o','markersize',mfs-3.4,'linewidth',lfs-1.5);
plot(ctsbma(~pptsbma(:,a)&pptsbma(:,b),a), ctsbma(~pptsbma(:,a)&pptsbma(:,b),b),'o','markersize',mfs-3.6,'linewidth',lfs-1.5);
title(sprintf('corrs time vs %s shuff %s',yl,t));
if a==2;plot(0,0, 'x');end %FOR PRE
plot(ctsbma(~pptsbma(:,a)&~pptsbma(:,b),a), ctsbma(~pptsbma(:,a)&~pptsbma(:,b),b),'o','markersize',mfs-3.8,'linewidth',lfs-1.5);
axis square; hold off;
xlabel('time'); ylabel(yl); set(gca, 'YLim',[-0.02,0.14]); set(gca, 'XLim',[-0.08,0.08]);
legend({'sig sig'; 'sig non'; 'non sig'; 'non non'},'location','best')

ctsbma(43,:)=[]; pptsbma(43,:)=[]; pairs(43,:)=[];
%}

%% A : Before During Temporal
gA = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.5,'A','fontweight','bold','fontsize',fs);
axA = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
plot(axA,ctsbma(:,1),            ctsbma(:,2),'.','markersize',mfs,'linewidth',lfs); hold on;
plot(axA,ctsbma(pptsbma(:,1),1), ctsbma(pptsbma(:,1),2),'go','markersize',mfs-3,'linewidth',lfs-1.5);
plot(axA,ctsbma(pptsbma(:,2),1), ctsbma(pptsbma(:,2),2),'mo','markersize',mfs,'linewidth',lfs-1.5);
f1 = fit(ctsbma(:,1),            ctsbma(:,2),'poly1');
mm = slimd(axA); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{}); 'a=%.2fx+%.2f'
%text(axA,0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'Units','normalized');
text(axA,0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(:,1),ctsbma(:,2)),2)),'Units','normalized');
xlabel('pre','fontsize',afs); 
ylabel('dur','fontsize',afs)
title('Time Correlations','fontsize',tfs); legend off;
axA.XLim = mm; axA.YLim = axA.XLim; axis(axA,'square');

set(gA,'Heights', [25 -1]);



%% B : Before During Spatial
gB = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.5,'B','fontweight','bold','fontsize',fs);
axB = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'),'visible','off');
plot(axB,ctsbma(:,4),            ctsbma(:,5),'.','markersize',mfs,'linewidth',lfs); hold on;
plot(axB,ctsbma(pptsbma(:,4),4), ctsbma(pptsbma(:,4),5),'go','markersize',mfs-3,'linewidth',lfs-1.5);
plot(axB,ctsbma(pptsbma(:,5),4), ctsbma(pptsbma(:,5),5),'mo','markersize',mfs,'linewidth',lfs-1.5);
f1 = fit(ctsbma(:,4),            ctsbma(:,5),'poly1');
mm = slimd(axB); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{});
%text(axB,0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'Units','normalized');
text(axB,0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(:,4),ctsbma(:,5)),2)),'Units','normalized');
xlabel('pre','fontsize',afs); 
ylabel('dur','fontsize',afs)
title('Spatial Correlations','fontsize',tfs); legend off;
axB.XLim = mm; axB.YLim = axB.XLim; axis(axB,'square');

set(gB,'Heights', [25 -1]);

%% BB : Before During XGS
gBB = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gBB,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.5,'C','fontweight','bold','fontsize',fs);
axBB = axes('Parent',uicontainer('Parent',gBB,'BackgroundColor','w'),'visible','off');
plot(axBB,ctsbma(:,7),            ctsbma(:,8),'.','markersize',mfs,'linewidth',lfs); hold on;
plot(axBB,ctsbma(pptsbma(:,7),7), ctsbma(pptsbma(:,7),8),'go','markersize',mfs-3,'linewidth',lfs-1.5);
plot(axBB,ctsbma(pptsbma(:,8),7), ctsbma(pptsbma(:,8),8),'mo','markersize',mfs,'linewidth',lfs-1.5);
f1 = fit(ctsbma(:,7),            ctsbma(:,8),'poly1');
mm = slimd(axBB); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{});
%text(axBB,0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'Units','normalized');
text(axBB,0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(:,7),ctsbma(:,8)),2)),'Units','normalized');
xlabel('pre','fontsize',afs); 
ylabel('dur','fontsize',afs)
title('XGS Correlations','fontsize',tfs); legend off;
axBB.XLim = mm; axBB.YLim = axBB.XLim; axis(axBB,'square');

set(gBB,'Heights', [25 -1]);


%% C : Before After Temporal
gC = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gC,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.5,'D','fontweight','bold','fontsize',fs);
axC = axes('Parent',uicontainer('Parent',gC,'BackgroundColor','w'),'visible','off');
plot(axC,ctsbma(:,1),            ctsbma(:,3),'.','markersize',mfs,'linewidth',lfs); hold on;
plot(axC,ctsbma(pptsbma(:,1),1), ctsbma(pptsbma(:,1),3),'go','markersize',mfs-3,'linewidth',lfs-1.5);
plot(axC,ctsbma(pptsbma(:,3),1), ctsbma(pptsbma(:,3),3),'ro','markersize',mfs,'linewidth',lfs-1.5);
f1 = fit(ctsbma(:,1),            ctsbma(:,3),'poly1');
mm = slimd(axC); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{});
%text(axC,0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'Units','normalized');
text(axC,0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(:,1),ctsbma(:,3)),2)),'Units','normalized');
xlabel('pre','fontsize',afs); 
ylabel('post','fontsize',afs)
title('Time Correlations','fontsize',tfs); legend off;
axC.XLim = mm; axC.YLim = axC.XLim; axis(axC,'square');

set(gC,'Heights', [25 -1]);


%% D : Before After Spatial
gD = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gD,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.5,'E','fontweight','bold','fontsize',fs);
axD = axes('Parent',uicontainer('Parent',gD,'BackgroundColor','w'),'visible','off');
plot(axD,ctsbma(:,4),            ctsbma(:,6),'.','markersize',mfs,'linewidth',lfs); hold on;
plot(axD,ctsbma(pptsbma(:,4),4), ctsbma(pptsbma(:,4),6),'go','markersize',mfs-3,'linewidth',lfs-1.5);
plot(axD,ctsbma(pptsbma(:,6),4), ctsbma(pptsbma(:,6),6),'ro','markersize',mfs,'linewidth',lfs-1.5);
f1 = fit(ctsbma(:,4),            ctsbma(:,6),'poly1');
mm = slimd(axD); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{});
%text(axD,0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'Units','normalized');
text(axD,0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(:,4),ctsbma(:,6)),2)),'Units','normalized');
xlabel('pre','fontsize',afs); 
ylabel('post','fontsize',afs)
title('Spatial Correlations','fontsize',tfs); legend off;
axD.XLim = mm; axD.YLim = axD.XLim; axis(axD,'square');
set(gD,'Heights', [25 -1]);

%% DD : Before After XGS
gDD = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gDD,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.5,'F','fontweight','bold','fontsize',fs);
axDD = axes('Parent',uicontainer('Parent',gDD,'BackgroundColor','w'),'visible','off');
plot(axDD,ctsbma(:,7),            ctsbma(:,9),'.','markersize',mfs,'linewidth',lfs); hold on;
plot(axDD,ctsbma(pptsbma(:,7),7), ctsbma(pptsbma(:,7),9),'go','markersize',mfs-3,'linewidth',lfs-1.5);
plot(axDD,ctsbma(pptsbma(:,9),7), ctsbma(pptsbma(:,9),9),'ro','markersize',mfs,'linewidth',lfs-1.5);
f1 = fit(ctsbma(:,7),            ctsbma(:,9),'poly1');
mm = slimd(axDD); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{});
%text(axDD,0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'Units','normalized');
text(axDD,0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(:,7),ctsbma(:,9)),2)),'Units','normalized');
xlabel('pre','fontsize',afs); 
ylabel('post','fontsize',afs)
title('XGS Correlations','fontsize',tfs); legend off;
axDD.XLim = mm; axDD.YLim = axDD.XLim; axis(axDD,'square');
set(gDD,'Heights', [25 -1])

% HISTOGRAMS
gE = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gE,'BackgroundColor','w'),'visible','off');
text(0,0.5,'G','fontweight','bold','fontsize',fs);
axE = axes('Parent',uicontainer('Parent',gE,'BackgroundColor','w'),'visible','off');
bar(axE,[sum(pptsbma(:,1))/len(pairs),0,0],'g');hold on;
bar(axE,[0,sum(pptsbma(:,2))/len(pairs),0],'m');
bar(axE,[0,0,sum(pptsbma(:,3))/len(pairs)],'r');
axE.YLim=[0,1];
legend({'pre';'dur';'post'}); 
title('Significant Time Correlations');
axis square; 
set(gE,'Heights', [25 -1]);

gF = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gF,'BackgroundColor','w'),'visible','off');
text(0,0.5,'H','fontweight','bold','fontsize',fs);
axF = axes('Parent',uicontainer('Parent',gF,'BackgroundColor','w'),'visible','off');
bar(axF,[sum(pptsbma(:,4))/len(pairs),0,0],'g');hold on;
bar(axF,[0,sum(pptsbma(:,5))/len(pairs),0],'m');
bar(axF,[0,0,sum(pptsbma(:,6))/len(pairs)],'r');
axF.YLim=[0,1];
legend({'pre';'dur';'post'}); 
title('Significant Spatial Correlations');
axis square; 
set(gF,'Heights', [25 -1]);

gFF = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gFF,'BackgroundColor','w'),'visible','off');
text(0,0.5,'I','fontweight','bold','fontsize',fs);
axFF = axes('Parent',uicontainer('Parent',gFF,'BackgroundColor','w'),'visible','off');
bar(axFF,[sum(pptsbma(:,7))/len(pairs),0,0],'g');hold on;
bar(axFF,[0,sum(pptsbma(:,8))/len(pairs),0],'m');
bar(axFF,[0,0,sum(pptsbma(:,9))/len(pairs)],'r');
axFF.YLim=[0,1];
legend({'pre';'dur';'post'}); 
title('Significant XGS Correlations');
axis square; 
set(gFF,'Heights', [25 -1]);

%% epilogue
set(gtop,'Widths', [-1 -1 -1]);%,'Widths', [-1 -1 -1]);

