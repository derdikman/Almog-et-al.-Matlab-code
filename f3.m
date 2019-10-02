%% preamble
% load('.\data\shuffling1000nanv2');

function f3(pairs,ctsbma, pptsbma)

dbstop if error  
nshuffle=1000;
fig = figure(993);clf; tic
set(fig,'Color','w', 'Position', [200 0 1000 650]); 
fss='fontsize'; mss='markersize';fw={'fontweight','bold'}; pnt='Parent';
fs ={fss 14}; afs ={fss 10}; mfs ={mss 10}; tfs = {fss,12};hts='Heights'; 


gTop = uix.GridFlex(pnt, fig,'Spacing',5,'padding',5, 'BackgroundColor','w','DividerMarkings','off');
gfp={'Spacing',5, 'BackgroundColor','w','DividerMarkings','off'};

tx=0;ty=0.5;tp={fw{:},fs{:}};
uip={'BackgroundColor','w'};
aup={'visible','off'};
%{
'sig:sig   sig:non   non:sig   non:non'
[sum(pptsbma(:,a)&pptsbma(:,b)),...
sum(pptsbma(:,a)&~pptsbma(:,b)),... 
sum(~pptsbma(:,a)&pptsbma(:,b)),...
sum(~pptsbma(:,a)&~pptsbma(:,b)) ]
if(a==1) disp('pre'); else disp('dur'); end
%}

%% A1 : Before During Temporal
gA1 = uix.GridFlex(pnt,gTop,gfp{:});
axes(pnt,uicontainer(pnt,gA1,uip{:}),aup{:});
text(tx,ty,'A',tp{:});
axes(pnt,uicontainer(pnt,gA1,uip{:}),aup{:});
x=ctsbma(:,1);y=ctsbma(:,2);arg=[];arg.show=[0 1 1];
plotARP(x,y,arg);
xlabel('pre',afs{:}); 
ylabel('dur',afs{:},fw{:});
title('time correlations',tfs{:}); legend off;

set(gA1,hts, [25 -1]);


%% C1 : Before During Spatial
gC1 = uix.GridFlex(pnt,gTop,gfp{:});
axes(pnt,uicontainer(pnt,gC1,uip{:}),aup{:});
text(tx,ty,'C',tp{:});
axes(pnt,uicontainer(pnt,gC1,uip{:}),aup{:});
x=ctsbma(:,4);y=ctsbma(:,5);arg=[];arg.show=[0 1 1];
plotARP(x,y,arg);
xlabel('pre',afs{:}); 
ylabel('dur',afs{:},fw{:});
title('spatial correlations',tfs{:}); legend off;

set(gC1,hts, [25 -1]);

%% A2 : Before After Temporal
gA2 = uix.GridFlex(pnt,gTop,gfp{:});
axes(pnt,uicontainer(pnt,gA2,uip{:}),aup{:});
text(tx,ty,'',tp{:});
axes(pnt,uicontainer(pnt,gA2,uip{:}),aup{:});
x=ctsbma(:,1);y=ctsbma(:,3);arg=[];arg.show=[0 1 1];
plotARP(x,y,arg);
xlabel('pre',afs{:}); 
ylabel('post',afs{:},fw{:});

set(gA2,hts, [25 -1]);


%% C2 : Before After Spatial
gC2 = uix.GridFlex(pnt,gTop,gfp{:});
axes(pnt,uicontainer(pnt,gC2,uip{:}),aup{:});
text(tx,ty,'',tp{:});
axes(pnt,uicontainer(pnt,gC2,uip{:}),aup{:});
x=ctsbma(:,4);y=ctsbma(:,6);arg=[];arg.show=[0 1 1];
plotARP(x,y,arg);
xlabel('pre',afs{:}); 
ylabel('post',afs{:},fw{:});

set(gC2,hts, [25 -1]);

N = nshuffle;
%% HISTOGRAMs
gB = uix.GridFlex(pnt,gTop,gfp{:});
axes(pnt,uicontainer(pnt,gB,uip{:}),aup{:});
text(0,0.0,'B',tp{:});
gBB = uix.GridFlex(pnt,gB,gfp{:});
axes(pnt,uicontainer(pnt,gBB,uip{:}),aup{:});
lp=len(pairs);
bar([sum(pptsbma(ctsbma(:,1)>0,1))/lp, sum(pptsbma(ctsbma(:,1)<0,1))/lp;...
     sum(pptsbma(ctsbma(:,2)>0,2))/lp, sum(pptsbma(ctsbma(:,2)<0,2))/lp;...
     sum(pptsbma(ctsbma(:,3)>0,3))/lp, sum(pptsbma(ctsbma(:,3)<0,3))/lp;]*100,'stacked')
ylabel('% total',fss,10);
set(gca,'xticklabel',{'pre';'dur';'post'},'YLim',[0,100],'xticklabelrotation',90);
title(sprintf('shuffling significance'),tfs{:});
legend({'(+) corr';'(-)  corr'}); axis square;
set(gB,hts, [25 -1]);

gD = uix.GridFlex(pnt,gTop,gfp{:});
axes(pnt,uicontainer(pnt,gD,uip{:}),aup{:});
text(0,0.0,'D',tp{:});
gDD = uix.GridFlex(pnt,gD,gfp{:});
axes(pnt,uicontainer(pnt,gDD,uip{:}),aup{:});
bar([sum(pptsbma(ctsbma(:,4)>0,4))/lp, sum(pptsbma(ctsbma(:,4)<0,4))/lp;...
     sum(pptsbma(ctsbma(:,5)>0,5))/lp, sum(pptsbma(ctsbma(:,5)<0,5))/lp;...
     sum(pptsbma(ctsbma(:,6)>0,6))/lp, sum(pptsbma(ctsbma(:,6)<0,6))/lp;]*100,'stacked')
ylabel('% total',fss,10);
set(gca,'xticklabel',{'pre';'dur';'post'},'YLim',[0,100],'xticklabelrotation',90);
title(sprintf('shuffling significance'),tfs{:});
legend({'(+) corr';'(-)  corr'}); axis square;

set(gD,hts, [25 -1]);
%space

%% epilogue
set(gTop,'Widths', [-1 -1 -1.5]);%,'Widths', [-1 -1 -1]);

clear gTop; clear gA1;clear gA2;clear gC1;clear gC2; clear gB; clear gD;   
clear fss; clear mss;clear fw;clear pnt;
clear fs;clear afs;clear mfs;clear tfs; clear nshuffle
clear tx;clear ty;clear tp;
clear uip;clear aup;
