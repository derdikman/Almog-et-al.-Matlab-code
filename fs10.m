%different grid threshes
% f3s5
function fs10(cellsn)
load('./data/thresh 3D nan');
%{
bt=0.9:-0.1:0.2;
mt=0.5:-0.1:0; 
tbm=[];sbm=[];tbmp=[];sbmp=[];tba=[];sba=[];tbap=[];sbap=[];
for b=bt
    v=[];w=[];o=[];u=[];vp=[];wp=[];op=[];up=[]; tic
    for m=mt
        [b m]
        prs=pairsMake(cellsn,b,m);
        t=tsbma(cellsn,prs);
        [r,p]=ccof(t(:,1),t(:,2));v(end+1)=r;vp(end+1)=p;
        [r,p]=ccof(t(:,4),t(:,5));w(end+1)=r;wp(end+1)=p;
        [r,p]=ccof(t(:,1),t(:,3));o(end+1)=r;op(end+1)=p;
        [r,p]=ccof(t(:,4),t(:,6));u(end+1)=r;up(end+1)=p;
    end
    tbm(end+1,:)=v;tbmp(end+1,:)=vp;
    sbm(end+1,:)=w;sbmp(end+1,:)=wp;
    tba(end+1,:)=o;tbap(end+1,:)=op;
    sba(end+1,:)=u;sbap(end+1,:)=up;
    toc
end
%}
x=repmat(bt',1,len(mt));y=repmat(mt,len(bt),1); 
p=(4-1)*len(bt)+5; ['thresh paper ' n2(x(p)) ' ' n2(y(p))]%0.5 0.2 thresh in paper
xl='min score pre'; yl='max score dur'; zl='corr';zlp=[0 0.2];

figure(9910); clf; set(gcf,'position',[10,10, 840,930]); 
tl=' [pre x dur]';
subplot(421);plot3(x,   y,   tbm,   'ko','markerfacecolor','g');
hold on;     plot3(x(p),y(p),tbm(p),'ko','markerfacecolor','b');
title(['temporal corrs by grid thresh'  tl]);grid on;
xlabel(xl);ylabel(yl);zlabel(zl);  zlim([0 1]);axis square;
subplot(422);plot3(x,   y,   tbmp,   'ko','markerfacecolor','g');
hold on;     plot3(x(p),y(p),tbmp(p),'ko','markerfacecolor','b');
grid on; title(['temporal corrs pval by grid thresh' tl]); 
xlabel(xl);ylabel(yl);zlabel('pval'); zlim(zlp);axis square;
subplot(425);plot3(x,   y,   sbm,   'ko','markerfacecolor','m'); 
hold on;     plot3(x(p),y(p),sbm(p),'ko','markerfacecolor','b');
grid on; title(['spatial corrs by grid thresh' tl]); 
xlabel(xl);ylabel(yl);zlabel(zl); zlim([0 1]);axis square;
subplot(426);plot3(x,   y,   sbmp,   'ko','markerfacecolor','m');
hold on;     plot3(x(p),y(p),sbmp(p),'ko','markerfacecolor','b');
grid on; title(['spatial corrs pval by grid thresh' tl]); 
xlabel(xl);ylabel(yl);zlabel('pval'); zlim(zlp);axis square;
%figure(4); clf; set(gcf,'position',[100,100, 900,900]);
tl=' [pre x post]';
subplot(423);plot3(x,   y,   tba,   'ko','markerfacecolor','c');
hold on;     plot3(x(p),y(p),tba(p),'ko','markerfacecolor','b');
title(['temporal corrs by grid thresh' tl]);grid on; zlim([0 1]); 
xlabel(xl);ylabel(yl);zlabel(zl); axis square;
subplot(424);plot3(x,   y,   tbap,   'ko','markerfacecolor','c');
hold on;     plot3(x(p),y(p),tbap(p),'ko','markerfacecolor','b');
grid on; title(['temporal corrs pval by grid thresh' tl]); 
xlabel(xl);ylabel(yl);zlabel('pval'); zlim(zlp);axis square;
subplot(427);plot3(x,   y,   sba,   'ko','markerfacecolor','r'); 
hold on;     plot3(x(p),y(p),sba(p),'ko','markerfacecolor','b');
grid on; title(['spatial corrs by grid thresh' tl]); 
xlabel(xl);ylabel(yl);zlabel(zl); zlim([0 1]);axis square;
subplot(428);plot3(x,   y,   sbap,   'ko','markerfacecolor','r');
hold on;     plot3(x(p),y(p),sbap(p),'ko','markerfacecolor','b');
grid on; title(['spatial corrs pval by grid thresh' tl]); 
xlabel(xl);ylabel(yl);zlabel('pval'); zlim(zlp);axis square;
suptitle('cell pair correlations and p-values across sessions by grid score threshold')
% figure(5);  clf; colormap default;
% subplot(121);
% surf(x,y,tbm);zlim([0 1]);hold on;  
% surf(x,y,tba);zlim([0 1]);
% subplot(122);
% surf(x,y,sbm);zlim([0 1]);hold on;  
% surf(x,y,sba);zlim([0 1]);

end