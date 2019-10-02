%{
% Load before calling function!!
load('C:\Noam\Data\muscimol\cells15nan');
load('.\data\shuffling1000nanv2');
load('.\data\pairsv2'); cels = unique(pairs(:))';
%}

function f4(cellsn, pairs,ctsbma)
tic
%% preamble
dbstop if error 

fig = figure(994); clf;  met='g'; %<<<<<'CHANGE FOR PAPER'
set(fig,'color','w', 'Position', [70 -300 1100 1300]); %1000 x

fss='fontsize';mss='markersize';p='Parent';fwss='fontweight';bss='bold';vss='visible';noss='none';
bcss='BackgroundColor'; btss='BorderType';spss='spacing';dmss='DividerMarkings';ofss='off';onss='on';
hzss='HorizontalAlignment';cnss='center'; coss='color';unss='Units';nzss='normalized'; htss='Heights';
s = {'before';'midall';'after'};tss = {'pre','dur','post'}; wtss='Widths';
vs='rayleigh_score';va='rayleigh_angle';
%fs ={fss 14}; afs ={fss 10}; mfs ={mss 10}; tfs = {fss,12}; fw={fwss,bss};
tx=0;ty=0.5; ty2=0;
fs = 12; afs = 10; mfs = 10; tfs = 12;  lfs = 3; h = 25; sp=1;

rst=0.4; rsth=rst; rstl=rst; %SAME AS PAIRS
ci = [112 114 117 118 119];%, 117 118 119]; %ci = [103, 104, 105, 100,111];  ci = [282,283,287]; ci = [112, 114, 115, 116, 117,118,119]; 115
cels = unique(pairs(:))';
clls = cellsn(cels); %t = [clls.midall]; t = [t.rayleigh_score]; hdclls = clls(t>rst); t=[];


%% make clusters
rsp=[];rs=[]; rmi=true(len(ctsbma),1);
for i=1:length(pairs)
    rsp(i,:) = [cellsn(pairs(i,1)).before.rayleigh_score,...
                cellsn(pairs(i,2)).before.rayleigh_score,...
                cellsn(pairs(i,1)).midall.rayleigh_score,...
                cellsn(pairs(i,2)).midall.rayleigh_score];
end
for i=1:length(cels)
    rs(i,:)  = [cellsn(cels(i)).before.rayleigh_score,...
        cellsn(cels(i)).midall.rayleigh_score];
end
%HD CELLS
hdc = rs(:,2)>=rsth & rs(:,1)<rstl;
%NON HD CELLS
hdcn = rs(:,2) <rstl & rs(:,1)<rstl; %ccl1 = ~ccl2; 
%rsp: rscore[c1b c2b c1m c2m]
%HD PAIRS
hdp = rsp(:,3)>=rsth & rsp(:,4)>=rsth ...   %midall
    & rsp(:,1)< rstl & rsp(:,2)< rstl;      %before
%NON HD PAIRS
hdpn = rsp(:,3)< rstl & rsp(:,4)< rstl ...     %midall
    & rsp(:,1)< rstl & rsp(:,2)< rstl;     %before
clear rsp; clear i; clear rsth; clear rstl; %clear rs;
chd = cels(hdc);
hdclls=cellsn(chd);
[sum(hdp) sum(hdpn)];
%rmi = ctsbma(:,1)<0.15 & ctsbma(:,2)<0.07; LOOK INTO INCLUDE IN PAPER

%% GRIDS TOP
gAll  = uix.GridFlex(p,  fig,spss,sp,bcss,'w',dmss,ofss);%,'padding',1);
gTop  = uix.GridFlex(p, gAll,spss,sp,bcss,'w',dmss,ofss);
gBot  = uix.GridFlex(p, gAll,spss,sp,bcss,'w',dmss,ofss);
gLeft = uix.GridFlex(p, gTop,spss,sp,bcss,'w',dmss,ofss);
gCent = uix.GridFlex(p, gTop,spss,sp,bcss,'w',dmss,ofss);
gRigt = uix.GridFlex(p, gTop,spss,sp,bcss,'w',dmss,ofss);



%% A
%axes(p,uicontainer(p,gleft,bcss,'w'),vss,ofss);

gA = uix.GridFlex(p,gLeft,spss,1,bcss,'w',dmss,onss);
axes(p,uix.Panel(p,gA,bcss,'w',btss,noss),vss,ofss);
text(tx,ty,'A',fwss,bss,fss,fs); axis off;
%%{
axA = {}; pf = []; ps = [];
for i = 1:len(ci)
    if i>1
       uix.Empty(p,gA);
    end  
    c = cellsn(ci(i));
    %TRAJ
    axes(p,uicontainer(p,gA,bcss,'w'),vss,ofss);
    text(0.5,0.5,sprintf('C%d',i),fwss,bss,fss,mfs,hzss,cnss);
            ps.x=' '; ps.y = ps.x; ps.t ='';
    axA{end+1} = axes(p,uicontainer(p,gA,bcss,'w'));
    axA{end}   = plotTR(axA{end},c.before,ps); ps.t = '';
    axA{end+1} = axes(p,uicontainer(p,gA,bcss,'w'));
    axA{end}   = plotTR(axA{end},c.midall,ps);
    %AC
            ps.t = sprintf('grid=%.1f',round(c.before.gridscore,2));
    axA{end+1} = axes(p,uicontainer(p,gA,bcss,'w'));
    axA{end}   = plotAC(axA{end},c.before,ps);
            ps.t = sprintf('grid=%.1f',round(c.midall.gridscore,2));
    axA{end+1} = axes(p,uicontainer(p,gA,bcss,'w'));
    axA{end}   = plotAC(axA{end},c.midall,ps);
    %RAY
            ps.t = sprintf('rayleigh=%.1f',round(c.before.rayleigh_score,2)); 
    axA{end+1} = axes(p,uicontainer(p,gA,bcss,'w'));
    axA{end}   = plotHD(axA{end},c.before,ps,1); axis square;
            ps.t = sprintf('rayleigh=%.1f',round(c.midall.rayleigh_score,2));
    axA{end+1} = axes(p,uicontainer(p,gA,bcss,'w'));
    axA{end}   = plotHD(axA{end},c.midall,ps,1); axis square;
end

set(gA,htss, [h 15 -1 -1 -1 -1 -1 -1]);%,wtss, [-1 -1 -1]);
%}

%% B


gB = uix.GridFlex(p,gBot,spss,1,bcss,'w',dmss,ofss);
%uix.Empty(p,gB); %EMPTY
axes(p,uix.Panel(p,gB,bcss,'w',btss,noss),vss,ofss);
text(0,ty,'B',fwss,bss,fss,fs,hzss,cnss);  axis off;
%%{
axB = {}; pf = []; ps = [];
ps.x=''; ps.y = ''; ps.t = ''; pf.lag = 1000; pf.movmean = 100;%pf.off = 0;
for i = 1:len(ci)-1
    for ii = i+1:len(ci)
        c1 = cellsn(ci(i));c2 = cellsn(ci(ii));
        if i+ii>3
            %uix.Empty(p,gB); %EMPTY
            uix.Empty(p,gB);
        end        
        %s1.t = sprintf('C%dxC%d',i,ii);
        
        axes(p,uicontainer(p,gB,bcss,'w'),vss,ofss);
        text(0.3,0.5,sprintf('C%dxC%d',i,ii),fwss,bss,fss,mfs,hzss,cnss);
        axB{end+1} = axes(p,uicontainer(p,gB,bcss,'w'));
        axB{end}  = plotTimeCorr(axB{end},c1.before,c2.before,pf,ps);
        axB{end+1} = axes(p,uicontainer(p,gB,bcss,'w'));
        ps.t = '';
        axB{end}  = plotTimeCorr(axB{end},c1.midall,c2.midall,pf,ps); 
        axB{end+1} = axes(p,uicontainer(p,gB,bcss,'w')); 
        axB{end}  = plotSpaceCorr(axB{end},c1.before,c2.before,pf,ps,met);
        axB{end+1} = axes(p,uicontainer(p,gB,bcss,'w')); 
        axB{end}  = plotSpaceCorr(axB{end},c1.midall,c2.midall,pf,ps,met);
        %uix.Empty(p,gB); %EMPTY
        
        
    end
end
set(gB,htss, [h 15 -1 -1 -1 -1]);%,wtss, [-1 -1 -1]);
%}

%gCenter
%gCent =  uix.GridFlex(p,gRigt,spss,sp,bcss,'w',dmss,ofss);



%uix.Empty(p,gCent); %EMPTY;

% %LEGEND
% axes(p,uicontainer(p,gCent,bcss,'w'),vss,ofss);
% plot(1,1,'mo'); hold on %,mss,mfs
% plot(1,1,'b.',mss,mfs);
% plot(1,1,'r.',mss,mfs);
% plot(1,1,'w.',mss,mfs+15);box off; axis off; %'Location','south'
% legend({'cohort';'cluster 1';'cluster 2'},fss,afs, 'position',[0.3 0.4 0.4 0.2]);

%% C hd single
gC = uix.GridFlex(p,gCent,spss,1,bcss,'w',dmss,ofss);
axes(p,uicontainer(p,gC,bcss,'w'),vss,ofss);
text(0,ty,'C',fwss,bss,fss,fs);
axes(p,uicontainer(p,gC,bcss,'w'),vss,ofss);
plot(rs(:,1),rs(:,2),'mo'); hold on;  %,mss,mfs
%hold on; plot(rs(rs(:,2)>0.5,1),rs(rs(:,2)>0.5,2),'r.',mss,mfs);
plot(rs(hdcn,1),rs(hdcn,2),'b.',mss,mfs);
plot(rs(hdc,1),rs(hdc,2),'r.',mss,mfs);
%title('rayleigh score by cell');
%legend({'cohort';'cluster 1';'cluster 2'},fss,afs,'Location','east');
xticks([tx,ty,1]);xticklabels([tx,ty,1]);
set(gca,'YTick',[tx,ty,1]);set(gca,'yTickLabel',[tx,ty,1],fss,afs)
title('rayleigh score per cell'); 
xlabel('pre',fss,afs); ylim([0 1]); xlim([0 1]);
ylabel('dur',fss,afs); axis square; hold off;
set(gC,htss, [h -1]);
 legend({'coh';'cl1';'cl2'},fss,afs, 'position',[0.8 0.5 0.2 0.1]);


%% D pairs time by hd
gD = uix.GridFlex(p,gCent,spss,1,bcss,'w',dmss,ofss);
axes(p,uicontainer(p,gD,bcss,'w'),vss,ofss);
%Text
text(0,ty2,'D',fwss,bss,fss,fs);
axes(p,uicontainer(p,gD,bcss,'w'),vss,ofss);
thdpn=hdpn&rmi; thdp=hdp&rmi;x1= ctsbma(thdpn,1);y1=ctsbma(thdpn,2);x2= ctsbma(thdp,1);y2=ctsbma(thdp,2);
fpp(x1,y1,x2,y2);
title('temporal correlations',fss,tfs); legend off;
function fpp(x1,y1,x2,y2)
plot(x1, y1,'b.',mss,mfs,'linewidth',lfs); hold on; %MISSING SOME??
plot(x2, y2,'r.',mss,mfs); 
f1 = fit(x1,y1,'poly1');
mm = slimd(gca); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
f2 = fit(x2,y2,'poly1');
mm = slimd(gca); pscat = plot(mm,f2(mm),'r:'); pscat.LineWidth = lfs-1;
[~, ~, t] = ccof(x1,y1);aa=rnd(f1.p1,2);text(0.1,0.9,...
sprintf('%s',t),coss,'b',unss,nzss); %'a=%.2f %s',rnd(f1.p1,2),t)
%sprintf('%s',t),coss,'b',unss,nzss); %'a=%.2f %s',rnd(f1.p1,2),t)
[~, ~, t] = ccof(x2,y2);aa=rnd(f2.p1,2);text(0.1,0.8,...
sprintf('%s',t),coss,'r',unss,nzss); %'a=%.2f %s',rnd(f1.p1,2),t)
%sprintf('%s',t),coss,'r',unss,nzss); %'a=%.2f %s',rnd(f2.p1,2),t)
xlabel('pre',fss,afs); ylabel('dur',fss,afs);xlim(mm);ylim(mm); axis square;
end
set(gD,htss, [h -1]);

%% D2 pairs spatial by hd %fig D2
gE = uix.GridFlex(p,gCent,spss,1,bcss,'w',dmss,ofss);
axes(p,uicontainer(p,gE,bcss,'w'),vss,ofss);
text(0,ty2,' ',fwss,bss,fss,fs); %empty
axes(p,uicontainer(p,gE,bcss,'w'),vss,ofss);
%plot(ctsbma(rmi,4),   ctsbma(rmi,5),'kx',mss,mfs,'linewidth',lfs-.5);hold on;
x1= ctsbma(thdpn,4);y1=ctsbma(thdpn,5);x2= ctsbma(thdp,4);y2=ctsbma(thdp,5);
fpp(x1,y1,x2,y2);
title('spatial correlations',fss,tfs); legend off;
set(gE,htss, [h -1]);

%% E hd cell %fig E
gF = uix.GridFlex(p,gRigt,spss,1,bcss,'w',dmss,ofss);
axes(p,uicontainer(p,gF,bcss,'w'),vss,ofss);
text(tx,ty2,'E',fwss,bss,fss,fs); 
yl=12;
axes(p,uicontainer(p,gF,bcss,'w'),vss,ofss); 
% title(sprintf('r-angle hd cells r-score > %0.1f',rst));
x = 2; t = []; a=fp(x, hdclls, rst,yl,'hd cluster dur','r'); len(t);
set(gF,htss, [h -1]);   
        
%% F hd cell by session %fig F
gG = uix.GridFlex(p,gRigt,spss,1,bcss,'w',dmss,ofss);
axes(p,uicontainer(p,gG,bcss,'w'),vss,ofss);
text(tx,ty2,'F',fwss,bss,fss,fs);
t1 = 'r-angle '; t2=' all';
yl=41;
axes(p,uicontainer(p,gG,bcss,'w'),vss,ofss); %PRE ALL
x = 1; t = [];fp(x,cellsn,rst,yl,[t1 tss{x} t2]); len(t); 
uicontainer(p,gG,bcss,'w'); %empty; %PRE DUR
axes(p,uicontainer(p,gG,bcss,'w'),vss,ofss);  
x = 2; t = [];fp(x,cellsn,rst,yl,[t1 tss{x} t2]); len(t);
%uicontainer(p,gG,bcss,'w'); %empty; %PRE POST
%axes(p,uicontainer(p,gG,bcss,'w'),vss,ofss);  
%x = 3; t = [];fp(x,cellsn,rst,yl,[t1 tss{x} t2]); len(t);        
function h=fp(x,cll,rst,yl,tit,c)
    cla;
    t = [cll.(s{x})]; tl=arrayfun(@(z) len(z.st),t); tl=tl>=100;   
    ta = ([t.(va)]); ts = [t.(vs)];  tis = ts>rst; t = ta(tis & tl);
    h=histogram(rad2deg(t),-180:30:180); xlim([-180 180]);  axis square; ylim([0 yl]);
    h.FaceColor='g'; if nargin==6; h.FaceColor=c; h.FaceAlpha=1; end 
    title(tit,fss,tfs); xlabel('r-angle'); ylabel(' ');
    l=[-180 180]; d=char(176); sf='%d%c';
    set(gca,'xtick',l,'xticklabel',{sprintf(sf,l(1),d); sprintf(sf,l(2),d)});
    text(0.1,0.95,sprintf('n=%d a=%.1f%c',len(t),meanang(t),176),unss,nzss);
end        
set(gG,htss, [h -1 h -1]);    


%% epilogue

set(gLeft,htss, [-1]);
set(gCent,wtss, [-1]);
set(gRigt,htss, [-1 -2]);
set(gTop,wtss, [-2 -1 -1]);
set(gBot,htss, [-1]);
set(gAll,htss, [-7 -5]);


toc

end
%

% plot(ctsbma(tl1,4), ctsbma(tl1,5),'b.',mss,mfs,'linewidth',lfs);hold on;
% plot(ctsbma(tl2,4), ctsbma(tl2,5),'r.',mss,mfs);
% f1 = fit(ctsbma(tl1,4),ctsbma(tl1,5),'poly1');
% mm = slimd(gca); pscat = plot(mm,f1(mm),'b:'); pscat.LineWidth = lfs-1;
% f2 = fit(ctsbma(tl2,4),ctsbma(tl2,5),'poly1');
% mm = slimd(gca); pscat = plot(mm,f2(mm),'r:'); pscat.LineWidth = lfs-1;
% text(0.1,0.9,...
% sprintf('r=%.2f a=%.2f',round(ccof(ctsbma(tl1,4),ctsbma(tl1,5)),2),round(f1.p1,2)),coss,'b',unss,nzss);
% text(0.1,0.8,...
% sprintf('r=%.2f a=%.2f',round(ccof(ctsbma(tl2,4),ctsbma(tl2,5)),2),round(f2.p1,2)),coss,'r',unss,nzss);
% xlabel('pre',fss,afs);ylabel('dur',fss,afs);xlim(mm); ylim(mm); axis square;