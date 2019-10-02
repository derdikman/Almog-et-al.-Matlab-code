%% preamble
dbstop if error  
%load('.\\data\\figdata');
fig = figure(993); tic
set(fig,'Color','w', 'Position', [200 0 1000 600]);
fs = 14; afs = 10; mfs = 10; tfs = 12; lfs = 3;  p = {};
gtop = uix.GridFlex('Parent', fig,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');

%{
a=1; b=4; %a=2; b=5;
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
a=1; b=4; yl = 'space'; t = 'PRE'; subplot(211)
a=2; b=5; yl = 'space'; t = 'DUR'; subplot(212)% %figure(99322);
%  plot(ctsbma(:,a),            ctsbma(:,b),'.','markersize',mfs,'linewidth',lfs); hold on;
plot(ctsbma(pptsbma(:,a)&pptsbma(:,b),a), ctsbma(pptsbma(:,a)&pptsbma(:,b),b),'o','markersize',mfs-3.2,'linewidth',lfs-1.5); hold on;
plot(ctsbma(pptsbma(:,a)&~pptsbma(:,b),a), ctsbma(pptsbma(:,a)&~pptsbma(:,b),b),'o','markersize',mfs-3.4,'linewidth',lfs-1.5);
plot(ctsbma(~pptsbma(:,a)&pptsbma(:,b),a), ctsbma(~pptsbma(:,a)&pptsbma(:,b),b),'o','markersize',mfs-3.6,'linewidth',lfs-1.5);
title(sprintf('corrs time vs %s shuff %s',yl,t));
%if a==2;plot(0,0, 'x');end %FOR PRE
plot(ctsbma(~pptsbma(:,a)&~pptsbma(:,b),a), ctsbma(~pptsbma(:,a)&~pptsbma(:,b),b),'o','markersize',mfs-3.8,'linewidth',lfs-1.5);
text(0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(:,a),ctsbma(:,b)),2)),'Units','normalized');
axis square; hold off; 
xlabel('time'); ylabel(yl); axis 'tight'; %set(gca, 'YLim',slimd(gca)); set(gca, 'XLim',slimd(gca));
legend({'sig sig'; 'sig non'; 'non sig'; 'non non'},'location','best')


%}

%% A : Before During Temporal
gA = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.5,'A','fontweight','bold','fontsize',fs);
axA = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
d1=ctsbma(:,1);d2=ctsbma(:,2);
plot(axA,ctsbma(:,1),            ctsbma(:,2),'.','markersize',mfs,'linewidth',lfs); hold on;
%plot(axA,ctsbma(pptsbma(:,1),1), ctsbma(pptsbma(:,1),2),'go','markersize',mfs-3,'linewidth',lfs-1.5);
%plot(axA,ctsbma(pptsbma(:,2),1), ctsbma(pptsbma(:,2),2),'mo','markersize',mfs,'linewidth',lfs-1.5);
f1 = fit(ctsbma(:,1),            ctsbma(:,2),'poly1');
mm = slimd(gca); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1; %slimd!!!
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{}); 'a=%.2fx+%.2f'
%text(axA,0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'Units','normalized');
text(axA,0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(:,1),ctsbma(:,2)),2)),'Units','normalized');
xlabel('pre','fontsize',afs); 
ylabel('dur','fontsize',afs,'fontweight','bold');
title('time correlations','fontsize',tfs); legend off;
axA.XLim = mm; axA.YLim = axA.XLim; axis(axA,'square');

set(gA,'Heights', [25 -1]);



%% B : Before During Spatial
gB = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.5,'C','fontweight','bold','fontsize',fs);
axB = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'),'visible','off');
plot(axB,ctsbma(:,4),            ctsbma(:,5),'.','markersize',mfs,'linewidth',lfs); hold on;
%plot(axB,ctsbma(pptsbma(:,4),4), ctsbma(pptsbma(:,4),5),'go','markersize',mfs-3,'linewidth',lfs-1.5);
%plot(axB,ctsbma(pptsbma(:,5),4), ctsbma(pptsbma(:,5),5),'mo','markersize',mfs,'linewidth',lfs-1.5);
f1 = fit(ctsbma(:,4),            ctsbma(:,5),'poly1');
mm = slimd(axB); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{});
%text(axB,0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'Units','normalized');
%text(axD,0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(:,4),ctsbma(:,5)),2)),'Units','normalized');
text(axB,0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(:,4),ctsbma(:,5)),2)),'Units','normalized');
xlabel('pre','fontsize',afs); 
ylabel('dur','fontsize',afs,'fontweight','bold');
title('spatial correlations','fontsize',tfs); legend off;
axB.XLim = mm; axB.YLim = axB.XLim; axis(axB,'square');

set(gB,'Heights', [25 -1]);

%% C : Before After Temporal
gC = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gC,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.5,' ','fontweight','bold','fontsize',fs);
axC = axes('Parent',uicontainer('Parent',gC,'BackgroundColor','w'),'visible','off');
plot(axC,ctsbma(:,1),            ctsbma(:,3),'.','markersize',mfs,'linewidth',lfs); hold on;
%plot(axC,ctsbma(pptsbma(:,1),1), ctsbma(pptsbma(:,1),3),'go','markersize',mfs-3,'linewidth',lfs-1.5);
%plot(axC,ctsbma(pptsbma(:,3),1), ctsbma(pptsbma(:,3),3),'ro','markersize',mfs,'linewidth',lfs-1.5);
f1 = fit(ctsbma(:,1),            ctsbma(:,3),'poly1');
mm = slimd(axC); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{});
%text(axC,0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'Units','normalized');
text(axC,0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(:,1),ctsbma(:,3)),2)),'Units','normalized');
xlabel('pre','fontsize',afs); 
ylabel('post','fontsize',afs,'fontweight','bold');
%title('Time Correlations','fontsize',tfs); legend off;
axC.XLim = mm; axC.YLim = axC.XLim; axis(axC,'square');

set(gC,'Heights', [25 -1]);


%% D : Before After Spatial
gD = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gD,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.5,' ','fontweight','bold','fontsize',fs);
axD = axes('Parent',uicontainer('Parent',gD,'BackgroundColor','w'),'visible','off');
plot(axD,ctsbma(:,4),            ctsbma(:,6),'.','markersize',mfs,'linewidth',lfs); hold on;
%plot(axD,ctsbma(pptsbma(:,4),4), ctsbma(pptsbma(:,4),6),'go','markersize',mfs-3,'linewidth',lfs-1.5);
%plot(axD,ctsbma(pptsbma(:,6),4), ctsbma(pptsbma(:,6),6),'ro','markersize',mfs,'linewidth',lfs-1.5);
f1 = fit(ctsbma(:,4),            ctsbma(:,6),'poly1');
mm = slimd(axD); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{});
%text(axD,0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'Units','normalized');
text(axD,0.1,0.9,sprintf('r=%.2f',round(ccof(ctsbma(:,4),ctsbma(:,6)),2)),'Units','normalized');
xlabel('pre','fontsize',afs); 
ylabel('post','fontsize',afs,'fontweight','bold');
%title('Spatial Correlations','fontsize',tfs); legend off;
axD.XLim = mm; axD.YLim = axD.XLim; axis(axD,'square');

set(gD,'Heights', [25 -1]);

%% HISTOGRAMs
gET = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gET,'BackgroundColor','w'),'visible','off');
%axes('Parent',uix.Panel('Parent',gE,'BackgroundColor','w','BorderType','none'),'visible','off');
text(0,0.0,'B','fontweight','bold','fontsize',fs,'HorizontalAlignment','center');
%text(0,0.5,'E','fontweight','bold','fontsize',fs);
gE = uix.GridFlex('Parent',gET,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gE,'BackgroundColor','w'),'visible','off');
%axes('Parent',uix.Panel('Parent',gE,'BackgroundColor','w','BorderType','none'),'visible','off');
n = max(ptsbma(:,1));
bar([sum(ptsbma(:,1)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(:,1))/length(pairs);...
     sum(ptsbma(:,2)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(:,2))/length(pairs);...
     sum(ptsbma(:,3)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(:,3))/length(pairs);]*100,'stacked')
ylabel('% total','fontsize',10);
set(gca,'xticklabel',{'pre';'during';'post'},'YLim',[0,100],'xticklabelrotation',90);
title(sprintf('shuffling significance'));
legend({'(-)  corr';'(+) corr'}); axis square;
%{
axes('Parent',uicontainer('Parent',gE,'BackgroundColor','w'),'visible','off');
n = max(zptsbma(:,1));
bar([sum(zptsbma(:,1)<= n*0.01)/length(pairs), sum(n*0.99 <= zptsbma(:,1))/length(pairs);...
     sum(zptsbma(:,2)<= n*0.01)/length(pairs), sum(n*0.99 <= zptsbma(:,2))/length(pairs);...
     sum(zptsbma(:,3)<= n*0.01)/length(pairs), sum(n*0.99 <= zptsbma(:,3))/length(pairs);]*100,'stacked')
ylabel(' ','fontsize',10);
set(gca,'xticklabel',{'pre';'during';'post'},'YLim',[0,100],'xticklabelrotation',90);
title(sprintf('random cell\nsignificant'));
%}
%suptitle('Temporal');
%bar([sum(pptsbma(:,1))/len(pairs),0,0],'g');hold on;
%bar([0,sum(pptsbma(:,2))/len(pairs),0],'m');
%bar([0,0,sum(pptsbma(:,3))/len(pairs)],'r');
%axE.YLim=[0,1];
%axis square;
%%%set(gE,'Widths', [-1 -1]);
set(gET,'Heights', [25 -1]);



gFT = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gFT,'BackgroundColor','w'),'visible','off');
%axes('Parent',uix.Panel('Parent',gE,'BackgroundColor','w','BorderType','none'),'visible','off');
text(0,0.0,'D','fontweight','bold','fontsize',fs,'HorizontalAlignment','center');
%text(0,0.5,'E','fontweight','bold','fontsize',fs);
gF = uix.GridFlex('Parent',gFT,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gF,'BackgroundColor','w'),'visible','off');
%axes('Parent',uix.Panel('Parent',gE,'BackgroundColor','w','BorderType','none'),'visible','off');
n = max(ptsbma(:,1));
bar([sum(ptsbma(:,4)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(:,4))/length(pairs);...
     sum(ptsbma(:,5)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(:,5))/length(pairs);...
     sum(ptsbma(:,6)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(:,6))/length(pairs);]*100,'stacked')
ylabel('% total','fontsize',10);
set(gca,'xticklabel',{'pre';'during';'post'},'YLim',[0,100],'xticklabelrotation',90);
title(sprintf('shuffling significance'));
legend({'(-)  corr';'(+) corr'}); axis square;
%{
axes('Parent',uicontainer('Parent',gF,'BackgroundColor','w'),'visible','off');
n = max(zptsbma(:,1));
bar([sum(zptsbma(:,4)<= n*0.01)/length(pairs), sum(n*0.99 <= zptsbma(:,4))/length(pairs);...
     sum(zptsbma(:,5)<= n*0.01)/length(pairs), sum(n*0.99 <= zptsbma(:,5))/length(pairs);...
     sum(zptsbma(:,6)<= n*0.01)/length(pairs), sum(n*0.99 <= zptsbma(:,6))/length(pairs);]*100,'stacked')
ylabel(' ','fontsize',10);
set(gca,'xticklabel',{'pre';'during';'post'},'YLim',[0,100],'xticklabelrotation',90);
title(sprintf('random cell\nsignificant'));
%}
%suptitle('Temporal');
%bar([sum(pptsbma(:,1))/len(pairs),0,0],'g');hold on;
%bar([0,sum(pptsbma(:,2))/len(pairs),0],'m');
%bar([0,0,sum(pptsbma(:,3))/len(pairs)],'r');
%axE.YLim=[0,1];
%axis square;
%%%set(gF,'Widths', [-1 -1]);
set(gFT,'Heights', [25 -1]);
%space

%% epilogue
set(gtop,'Widths', [-1 -1 -1.5]);%,'Widths', [-1 -1 -1]);


%{
gF = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gF,'BackgroundColor','w'),'visible','off');
text(0,0.5,'F','fontweight','bold','fontsize',fs);
axF = axes('Parent',uicontainer('Parent',gF,'BackgroundColor','w'),'visible','off');
bar([sum(pptsbma(:,4))/len(pairs),0,0],'g');hold on;
bar([0,sum(pptsbma(:,5))/len(pairs),0],'m');
bar([0,0,sum(pptsbma(:,6))/len(pairs)],'r');
axF.YLim=[0,1];
legend({'pre';'dur';'post'}); 
title('Significant Space Correlations');
axis square; 
set(gF,'Heights', [25 -1]);


%}
