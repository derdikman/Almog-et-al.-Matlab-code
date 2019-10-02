% uncomment to load cells
%load(sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45));
%load('Z:\\data\\noam\\Data\\muscimol\\noam\\cells_45min_d_patchtraj_rayleigh'); %lab
load('.\data\pairs'); cels = unique(pairs(:))';
load('.\data\cptsbma'); rmi=[];

%% preamble
dbstop if error
%%loads();
%%pairs(43,:)=[];
%%{
rsp=[];rs = [];
for i=1:length(pairs)
    rsp(i,:) = [cells{pairs(i,1)}.before.rayleigh_score,...
        cells{pairs(i,2)}.before.rayleigh_score,...
        cells{pairs(i,1)}.midall.rayleigh_score,...
        cells{pairs(i,2)}.midall.rayleigh_score];
end
for i=1:length(cels)
    rs(i,:)  = [cells{cels(i)}.before.rayleigh_score,...
        cells{cels(i)}.midall.rayleigh_score];
end
cl2 = rsp(:,3)>0.4 & rsp(:,4)>0.4; cl1 = ~cl2;

% pairs(43,:)=[];ctsbma(43,:)=[]; ptsbma(43,:)=[];pptsbma(43,:)=[];cl1(43)=[];cl2(43)=[];zpptsbma(43,:)=[];
%}
f = figure(994); tic
set(f,'Color','w', 'Position', [200 0 1300 700]);
fs = 12; afs = 10; mfs = 10; tfs = 12; lfs = 3; p = {};

%%
gtop = uix.GridFlex('Parent', f,'Spacing',1, 'BackgroundColor','w','DividerMarkings','off'); 
w=25;

%% A hd single
gA = uix.GridFlex('Parent',gtop,'Spacing',1, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.9,'A','fontweight','bold','fontsize',fs);
axA = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
plot(rs(:,1),rs(:,2),'b.','markersize',mfs);
hold on; plot(rs(rs(:,2)>0.5,1),rs(rs(:,2)>0.5,2),'r.','markersize',mfs);
%title('rayleigh score by cell');
legend({'cluster 1';'cluster 2'},'fontsize',afs,'Location','east');
set(gca,'XTick',[0,0.5,1]);set(gca,'XTickLabel',[0,0.5,1],'fontsize',afs)
set(gca,'YTick',[0,0.5,1]);set(gca,'yTickLabel',[0,0.5,1],'fontsize',afs)
title('rayleigh score per cell');
xlabel('pre','fontsize',afs); axA.YLim=[0 1];axA.XLim=[0 1];
ylabel('dur','fontsize',afs);axis(axA,'square');
set(gA,'Heights', [w -1]);

%% D AVERAGE
gFT = uix.GridFlex('Parent',gtop,'Spacing',1, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gFT,'BackgroundColor','w'),'visible','off');
text(0,0.9,'D','fontweight','bold','fontsize',fs);
axes('Parent',uicontainer('Parent',gFT,'BackgroundColor','w'),'visible','off');
text(0.5,0.0,'mean absolute correlation','fontsize',fs,'HorizontalAlignment','center');
gF = uix.GridFlex('Parent',gFT,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axA = axes('Parent',uicontainer('Parent',gF,'BackgroundColor','w'),'visible','off');
b=plot(axA,[mean(abs(ctsbma(:,1))), mean(abs(ctsbma(cl1,1))), mean(abs(ctsbma(cl2,1)));...
       mean(abs(ctsbma(:,2))), mean(abs(ctsbma(cl1,2))), mean(abs(ctsbma(cl2,2)));...
       mean(abs(ctsbma(:,3))), mean(abs(ctsbma(cl1,1))), mean(abs(ctsbma(cl2,3)))]);%,'grouped');
b(1).Color = 'k';b(2).Color = 'b';b(3).Color = 'r'; %FaceColor
set(axA,'xticklabel',{'pre';'during';'post'},'xticklabelrotation',90,'YLim',slimdy(axA));
ylabel('');title(sprintf('temporal'),'fontsize',fs-2); %axis square;
legend({'all';'cluster 1';'cluster 2'},'Location','north');
axA = axes('Parent',uicontainer('Parent',gF,'BackgroundColor','w'),'visible','off');
b=plot([mean(abs(ctsbma(:,4))), mean(abs(ctsbma(cl1,4))), mean(abs(ctsbma(cl2,4)));...
       mean(abs(ctsbma(:,5))), mean(abs(ctsbma(cl1,5))), mean(abs(ctsbma(cl2,5)));...
       mean(abs(ctsbma(:,6))), mean(abs(ctsbma(cl1,6))), mean(abs(ctsbma(cl2,6)))]);%,'grouped');
b(1).Color = 'k';b(2).Color = 'b';b(3).Color = 'r'; %FaceColor
set(gca,'xticklabel',{'pre';'during';'post'},'xticklabelrotation',90,'YLim',slimdy(gca));
title(sprintf('spatial'),'fontsize',fs-2); % ylabel('abs(correlation)'); %axis square;
%legend({'all';'cluster 1';'cluster 2'});
set(gF,'Widths', [-1 -1]);
set(gFT,'Heights', [25 25 -1]);


%% B time by hd
gB = uix.GridFlex('Parent',gtop,'Spacing',1, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.9,'B','fontweight','bold','fontsize',fs);
axB = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'),'visible','off');
%plot(ctsbma(rmi,1),   ctsbma(rmi,2),'kx','markersize',mfs,'linewidth',lfs-.5);hold on;  %REMOVE OUTLIER
plot(ctsbma(:,1),   ctsbma(:,2),'b.','markersize',mfs,'linewidth',lfs); hold on;
plot(ctsbma(cl2,1), ctsbma(cl2,2),'r.','markersize',mfs);
tc1=cl1;  tc2=cl2; %tc1(rmi)=false; tc2(rmi)=false; %REMOVE OUTLIER
f1 = fit(ctsbma(tc1,1),ctsbma(tc1,2),'poly1');
mm = slimd(gca); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
f2 = fit(ctsbma(tc2,1),ctsbma(tc2,2),'poly1');
mm = slimd(gca); pscat = plot(mm,f2(mm),'r:'); pscat.LineWidth = lfs-1;
%plot(ctsbma(pptsbma(:,2),1), ctsbma(pptsbma(:,2),2),'mo','markersize',mfs,'linewidth',lfs-2.5);
%plot(ctsbma(pptsbma(:,1),1), ctsbma(pptsbma(:,1),2),'go','markersize',mfs,'linewidth',lfs-2.5);
%f3 = fit(ctsbma(pptsbma(:,2),1), ctsbma(pptsbma(:,2),2),'poly1');
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{}); 'a=%.2fx+%.2f'
% text(0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'color','b','Units','normalized');
% text(0.1,0.8,sprintf('%.2fx + %.2f',round(f2.p1,2),round(f2.p2,2)),'color','r','Units','normalized');
% text(0.1,0.7,sprintf('%.2fx + %.2f',round(f3.p1,2),round(f3.p2,2)),'color','m','Units','normalized');
text(0.1,0.9,sprintf('r=%.2f a=%.2f',round(ccof(ctsbma(tc1,1),ctsbma(tc1,2)),2),round(f1.p1,2)),...
    'color','b','Units','normalized');
text(0.1,0.8,sprintf('r=%.2f a=%.2f',round(ccof(ctsbma(tc2,1),ctsbma(tc2,2)),2),round(f2.p1,2)),...
    'color','r','Units','normalized');
%text(0.1,0.7,sprintf('r=%.2f',round(ccof(ctsbma(pptsbma(:,2),1),ctsbma(pptsbma(:,2),2)),2)),'color','m','Units','normalized');
xlabel('pre','fontsize',afs);
ylabel('dur','fontsize',afs)
title('temporal correlations','fontsize',tfs); legend off;
axB.XLim = mm; axB.YLim = axB.XLim; axis(axB,'square');
set(gB,'Heights', [w -1]);

%% HISTOGRAM Significance
gDT = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gDT,'BackgroundColor','w'),'visible','off');
text(0,0.0,'E','fontweight','bold','fontsize',fs,'HorizontalAlignment','center');
axes('Parent',uicontainer('Parent',gDT,'BackgroundColor','w'),'visible','off');
text(0.5,0.0,'significant temporal correlations','fontsize',fs,'HorizontalAlignment','center');
%text(0.1,0.0,'significance(dur)','fontsize',fs-1.5,'HorizontalAlignment','center');
gD = uix.GridFlex('Parent',gDT,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gD,'BackgroundColor','w'),'visible','off');
n = max(ptsbma(:,1));
bar([sum(ptsbma(cl1,1)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(cl1,1))/length(pairs);...
     sum(ptsbma(cl1,2)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(cl1,2))/length(pairs);...
     sum(ptsbma(cl1,3)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(cl1,3))/length(pairs);]*100,'stacked')
ylabel('% total','fontsize',10);
set(gca,'xticklabel',{'pre';'during';'post'},'YLim',[0,100],'xticklabelrotation',90);
title(sprintf('cluster 1'),'color','b','fontsize',fs-2);
legend({'(-)  corr';'(+) corr'},'Location','north'); %axis square;
axes('Parent',uicontainer('Parent',gD,'BackgroundColor','w'),'visible','off');
n = max(ptsbma(:,1));
bar([sum(ptsbma(cl2,1)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(cl2,1))/length(pairs);...
     sum(ptsbma(cl2,2)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(cl2,2))/length(pairs);...
     sum(ptsbma(cl2,3)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(cl2,3))/length(pairs);]*100,'stacked')
ylabel(' ','fontsize',10);
set(gca,'xticklabel',{'pre';'during';'post'},'YLim',[0,100],'xticklabelrotation',90);
title(sprintf('cluster 2'),'color','r','fontsize',fs-2); %axis square;
set(gD,'Widths', [-1 -1]);
set(gDT,'Heights', [25 25 -1]);



%% C spatial by hd
gC = uix.GridFlex('Parent',gtop,'Spacing',1, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gC,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.9,'C','fontweight','bold','fontsize',fs);
axC = axes('Parent',uicontainer('Parent',gC,'BackgroundColor','w'),'visible','off');
%plot(ctsbma(rmi,4),   ctsbma(rmi,5),'kx','markersize',mfs,'linewidth',lfs-.5);hold on;
plot(ctsbma(:,4),   ctsbma(:,5),'b.','markersize',mfs,'linewidth',lfs);hold on;
plot(ctsbma(cl2,4), ctsbma(cl2,5),'r.','markersize',mfs);
tc1=cl1; tc1(rmi)=false; tc2=cl2; tc2(rmi)=false;
f1 = fit(ctsbma(cl1,4),ctsbma(cl1,5),'poly1');
mm = slimd(gca); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
f2 = fit(ctsbma(cl2,4),ctsbma(cl2,5),'poly1');
mm = slimd(gca); pscat = plot(mm,f2(mm),'r:'); pscat.LineWidth = lfs-1;
% plot(ctsbma(:,4),   ctsbma(:,5),'b.','markersize',mfs,'linewidth',lfs); hold on;
% plot(ctsbma(cl2,4), ctsbma(cl2,5),'r.','markersize',mfs);
% f1 = fit(ctsbma(cl1,4),ctsbma(cl1,5),'poly1');
% mm = slimd(gca); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
% f2 = fit(ctsbma(cl2,4),ctsbma(cl2,5),'poly1');
% mm = slimd(gca); pscat = plot(mm,f2(mm),'r:'); pscat.LineWidth = lfs-1;

%plot(ctsbma(pptsbma(:,5),4), ctsbma(pptsbma(:,5),5),'mo','markersize',mfs,'linewidth',lfs-2.5);
%f3 = fit(ctsbma(pptsbma(:,5),4), ctsbma(pptsbma(:,5),5),'poly1');
%textbest(axC,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),{},{}); 'a=%.2fx+%.2f'
% text(0.1,0.9,sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),'color','b','Units','normalized');
% text(0.1,0.8,sprintf('%.2fx + %.2f',round(f2.p1,2),round(f2.p2,2)),'color','r','Units','normalized');
% text(0.1,0.7,sprintf('%.2fx + %.2f',round(f3.p1,2),round(f3.p2,2)),'color','m','Units','normalized');
text(0.1,0.9,sprintf('r=%.2f a=%.2f',round(ccof(ctsbma(cl1,4),ctsbma(cl1,5)),2),round(f1.p1,2)),...
    'color','b','Units','normalized');
text(0.1,0.8,sprintf('r=%.2f a=%.2f',round(ccof(ctsbma(cl2,4),ctsbma(cl2,5)),2),round(f2.p1,2)),...
    'color','r','Units','normalized');
%text(0.1,0.7,sprintf('r=%.2f',round(ccof(ctsbma(pptsbma(:,5),4),ctsbma(pptsbma(:,5),5)),2)),'color','m','Units','normalized');
xlabel('pre','fontsize',afs);
ylabel('dur','fontsize',afs)
title('spatial correlations','fontsize',tfs); legend off;
axC.XLim = mm; axC.YLim = axC.XLim; axis(axC,'square');
set(gC,'Heights', [w -1]);


%% HISTOGRAMs
gET = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gET,'BackgroundColor','w'),'visible','off');
text(0,0.0,'F','fontweight','bold','fontsize',fs,'HorizontalAlignment','center');
axes('Parent',uicontainer('Parent',gET,'BackgroundColor','w'),'visible','off');
text(0.5,0.0,'significant spatial correlations','fontsize',fs,'HorizontalAlignment','center');

%text(0.5,0.0,'significance dur','fontsize',fs-1,'HorizontalAlignment','center');
gE = uix.GridFlex('Parent',gET,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axes('Parent',uicontainer('Parent',gE,'BackgroundColor','w'),'visible','off');
n = max(ptsbma(:,1));
bar([sum(ptsbma(cl1,4)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(cl1,4))/length(pairs);...
     sum(ptsbma(cl1,5)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(cl1,5))/length(pairs);...
     sum(ptsbma(cl1,6)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(cl1,6))/length(pairs);]*100,'stacked')
ylabel('% total','fontsize',10);
set(gca,'xticklabel',{'pre';'during';'post'},'YLim',[0,100],'xticklabelrotation',90);
title(sprintf('cluster 1'),'color','b','fontsize',fs-2); %axis square;
%legend({'(-)  corr';'(+) corr'}); 
axes('Parent',uicontainer('Parent',gE,'BackgroundColor','w'),'visible','off');
n = max(ptsbma(:,1));
bar([sum(ptsbma(cl2,1)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(cl2,1))/length(pairs);...
     sum(ptsbma(cl2,2)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(cl2,2))/length(pairs);...
     sum(ptsbma(cl2,3)<= n*0.01)/length(pairs), sum(n*0.99 <= ptsbma(cl2,3))/length(pairs);]*100,'stacked')
ylabel(' ','fontsize',10);
set(gca,'xticklabel',{'pre';'during';'post'},'YLim',[0,100],'xticklabelrotation',90);
title(sprintf('cluster 2'),'color','r','fontsize',fs-2);%axis square;
set(gE,'Widths', [-1 -1]);
set(gET,'Heights', [25 25 -1]);



%%
set(gtop,'Widths', [-1 -1 -1],'Heights', [-3 -2]);

