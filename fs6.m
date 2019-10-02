% grid score corr
% f3s2
%load('C:\Noam\Data\muscimol\cells15nan')
%load('.\data\shuffling1000nanv2')
%load('.\data\pairsv2')
%gridscore='gridscore'; 
%for v3%
%gridscore='gs2';
function fs6(cellsn,pairs,ctsbma)
%GRIDSCORE VS CORR
gridscore='gridscore';
gb=[];gd=[];
for j=1:len(pairs) % col1 min col2 max col3 mean
    c1 = cellsn(pairs(j,1)); c2 = cellsn(pairs(j,2));%c1.ind, c2.ind
    s1 = c1.before.(gridscore); s2 = c2.before.(gridscore); %s1(s1==-2)=0; s2(s2==-2)=0;
    gb(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
    s1 = c1.midall.(gridscore); s2 = c2.midall.(gridscore); %s1(s1==-2)=0; s2(s2==-2)=0;
    gd(j,:) = [min(s1,s2) max(s1,s2) mean([s1 s2])];
end


fig=figure(1006); clf; set(fig,'color','w', 'Position', [200 70 1100 550]);
arg=[];arg.setlim=false; xl='score'; yl='corr'; 

tos=1;%time pre       col3 mean  mb = before, max mean min of each pair
subplot(241);x = gb(:,3); y = ctsbma(:,tos); plotARP(x,y,arg);
title('time pre');xlabel(xl); ylabel(yl);
subplot(242);x = gd(:,3); y = ctsbma(:,tos+1);plotARP(x,y,arg);
title('time dur');xlabel(xl); ylabel(yl); 

tos = 4;%spatial pre
subplot(243);x = gb(:,3); y = ctsbma(:,tos); plotARP(x,y,arg);
title('spatial pre');xlabel(xl); ylabel(yl); 
subplot(244);x = gd(:,3); y = ctsbma(:,tos+1);plotARP(x,y,arg);
title('spatial dur');xlabel(xl); ylabel(yl); 

tos=1; yl = 'abs(corr)';%time pre
subplot(245);x = gb(:,3); y = ctsbma(:,tos); y = abs(y); plotARP(x,y,arg);
title('time pre');xlabel(xl); ylabel(yl); 
subplot(246);x = gd(:,3); y = ctsbma(:,tos+1); y = abs(y); plotARP(x,y,arg); 
title('time dur');xlabel(xl); ylabel(yl); 

tos = 4; %spatial pre
subplot(247);x = gb(:,3); y = ctsbma(:,tos); y = abs(y); plotARP(x,y,arg);
title('spatial pre');xlabel(xl); ylabel(yl); 
subplot(248);x = gd(:,3); y = ctsbma(:,tos+1); y = abs(y); plotARP(x,y,arg); 
title('spatial dur');xlabel(xl); ylabel(yl); 

suptitle('Mean grid score vs correlation');
end