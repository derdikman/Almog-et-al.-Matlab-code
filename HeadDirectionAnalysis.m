

trandclusters = [];
n=100000;tic
for i = 1:n
    rphd = logical(zeros(len(pairs),1));
    rnphd = logical(zeros(len(pairs),1));
    t=randperm(len(pairs),len(phd) + len(nphd))'; 
    rphd(t(1:len(phd)))=true;
    rnphd(t(len(phd)+1:end))=true;
    
    trandclusters(i,:)=[ccof(ctsbma(rphd,1),ctsbma(rphd,2)),...
                        ccof(ctsbma(rnphd,1),ctsbma(rnphd,2)),...
                        ccof(ctsbma(rphd,4),ctsbma(rphd,5)),...
                        ccof(ctsbma(rnphd,4),ctsbma(rnphd,5))];
end
toc
figure(1);
subplot(211)
stats(trandclusters(:,1)-trandclusters(:,2))
t = ccof(ctsbma(phd,1),ctsbma(phd,2)) - ccof(ctsbma(nphd,1),ctsbma(nphd,2))
histogram(trandclusters(:,1)-trandclusters(:,2),100)
line([t, t], ylim, 'LineWidth', 2, 'Color', 'r');
title('time')
subplot(212)
stats(trandclusters(:,3)-trandclusters(:,4))
t=ccof(ctsbma(phd,4),ctsbma(phd,5)) - ccof(ctsbma(nphd,4),ctsbma(nphd,5))
histogram(trandclusters(:,3)-trandclusters(:,4),100)
line([t, t], ylim, 'LineWidth', 2, 'Color', 'r');
title('spatial')
suptitle(sprintf('histogram of n=%d runs of random clusters \n[nonHD - HD] of corr(before,during)',n)) 


figure(1);
trandclusters = [];
n=10000
for i = 1:n
    rcl1 = logical(zeros(len(pairs),1));
    rcl1(randperm(len(pairs),sum(cl1)))=true;   
    trandclusters(i,:)=[ccof(ctsbma(rcl1,1),ctsbma(rcl1,2)),...
                        ccof(ctsbma(~rcl1,1),ctsbma(~rcl1,2)),...
                        ccof(ctsbma(rcl1,4),ctsbma(rcl1,5)),...
                        ccof(ctsbma(~rcl1,4),ctsbma(~rcl1,5))];
end

subplot(211)
%stats(trandclusters(:,1)-trandclusters(:,2))
t = ccof(ctsbma(cl1,1),ctsbma(cl1,2)) - ccof(ctsbma(~cl1,1),ctsbma(~cl1,2))
histogram(trandclusters(:,1)-trandclusters(:,2),20)
line([t, t], ylim, 'LineWidth', 2, 'Color', 'r');
title('time')

subplot(212)
%stats(trandclusters(:,3)-trandclusters(:,4))
t=ccof(ctsbma(cl1,4),ctsbma(cl1,5)) - ccof(ctsbma(~cl1,4),ctsbma(~cl1,5))
histogram(trandclusters(:,3)-trandclusters(:,4),20)
line([t, t], ylim, 'LineWidth', 2, 'Color', 'r');
title('spatial')

suptitle(sprintf('histogram of n=%d runs of random clusters \n[nonHD - HD] of corr(before,during)',n)) 

%%%%%%%%%%% SLOPE


crandclusters = [];
arandclusters = [];
nn=10000
for i = 1:nn
    rcl1 = logical(zeros(len(pairs),1));
    rcl1(randperm(len(pairs),sum(cl1)))=true;    
    crandclusters(i,:)=[ccof(ctsbma(rcl1,1),ctsbma(rcl1,2)),...
                        ccof(ctsbma(~rcl1,1),ctsbma(~rcl1,2)),...
                        ccof(ctsbma(rcl1,4),ctsbma(rcl1,5)),...
                        ccof(ctsbma(~rcl1,4),ctsbma(~rcl1,5))];
%      f1 = fit(ctsbma( rcl1,1),ctsbma( rcl1,2),'poly1');
%      f2 = fit(ctsbma(~rcl1,1),ctsbma(~rcl1,2),'poly1');
%      f3 = fit(ctsbma( rcl1,4),ctsbma( rcl1,5),'poly1');
%      f4 = fit(ctsbma(~rcl1,4),ctsbma(~rcl1,5),'poly1');
%      arandclusters(i,:) = [f1.p1, f2.p1, f3.p1, f4.p1];
     f1 = polyfit(ctsbma( rcl1,1),ctsbma( rcl1,2),1);
     f2 = polyfit(ctsbma(~rcl1,1),ctsbma(~rcl1,2),1);
     f3 = polyfit(ctsbma( rcl1,4),ctsbma( rcl1,5),1);
     f4 = polyfit(ctsbma(~rcl1,4),ctsbma(~rcl1,5),1);
     arandclusters(i,:) = [f1(1), f2(1), f3(1), f4(1)];
 i
end

figure(2);
subplot(221)
%stats(crandclusters(:,1)-crandclusters(:,2))
t = ccof(ctsbma(cl1,1),ctsbma(cl1,2)) - ccof(ctsbma(~cl1,1),ctsbma(~cl1,2));
tt = crandclusters(:,1)-crandclusters(:,2)
histogram(tt,20);
line([t, t], ylim, 'LineWidth', 2, 'Color', 'r');
sum(tt<t)/n
title(['time corr diff p=' num2str(sum(tt<t)/n)])

subplot(222)
%stats(crandclusters(:,3)-crandclusters(:,4))
t=ccof(ctsbma(cl1,4),ctsbma(cl1,5)) - ccof(ctsbma(~cl1,4),ctsbma(~cl1,5));
tt =crandclusters(:,3)-crandclusters(:,4);
histogram(tt,20);
line([t, t], ylim, 'LineWidth', 2, 'Color', 'r');
sum(tt<t)/n
title(['spatial corr diff p=' num2str(sum(tt<t)/n)])

subplot(223)
%stats(arandclusters(:,1)-arandclusters(:,2))
f1 = fit(ctsbma( cl1,1),ctsbma( cl1,2),'poly1');
f2 = fit(ctsbma(~cl1,1),ctsbma(~cl1,2),'poly1');
t = f1.p1 - f2.p1;
tt = arandclusters(:,1)-arandclusters(:,2);
histogram(tt,20);
line([t, t], ylim, 'LineWidth', 2, 'Color', 'r');
sum(tt<t)/n
title(['time slope diff p=' num2str(sum(tt<t)/n)])

subplot(224)
%stats(arandclusters(:,3)-arandclusters(:,4))
f3 = fit(ctsbma( cl1,4),ctsbma( cl1,5),'poly1');
f4 = fit(ctsbma(~cl1,4),ctsbma(~cl1,5),'poly1');
t = f3.p1 - f4.p1;
tt = arandclusters(:,3)-arandclusters(:,4);
histogram(tt,20);
line([t, t], ylim, 'LineWidth', 2, 'Color', 'r');
sum(tt<t)/n
title(['spatial slope diff p=' num2str(sum(tt<t)/n)])

suptitle(sprintf('histogram of n=%d runs of random clusters \n[before vs during] [nonHD - HD]',n)) 


%%%% ABS DIFF MEAN

trandclusters = [];
n=10000
for i = 1:n
    rcl1 = logical(zeros(len(pairs),1));
    rcl1(randperm(len(pairs),sum(cl1)))=true;
    
    trandclusters(i,:)=[mean(ctsbma(rcl1,1))-mean(ctsbma(rcl1,2)),...
                        mean(ctsbma(~rcl1,1))-mean(ctsbma(~rcl1,2)),...
                        mean(ctsbma(rcl1,4))-mean(ctsbma(rcl1,5)),...
                        mean(ctsbma(~rcl1,4))-mean(ctsbma(~rcl1,5))];
end

figure(3);
subplot(211)
%stats(trandclusters(:,1)-trandclusters(:,2))
t = mean(ctsbma(cl1,1))-mean(ctsbma(cl1,2)) - (mean(ctsbma(~cl1,1))-mean(ctsbma(~cl1,2)) );
tt = trandclusters(:,1)-trandclusters(:,2)
histogram(tt,100);
line([t, t], ylim, 'LineWidth', 2, 'Color', 'r');
sum(tt<t)/n
title(['time mean diff p=' num2str(sum(tt<t)/n)])

subplot(212)
%stats(trandclusters(:,3)-trandclusters(:,4))
t=mean(ctsbma(cl1,4))-mean(ctsbma(cl1,5)) - (mean(ctsbma(~cl1,4))-mean(ctsbma(~cl1,5)) );
tt =trandclusters(:,3)-trandclusters(:,4);
histogram(tt,100);
line([t, t], ylim, 'LineWidth', 2, 'Color', 'r');
sum(tt<t)/n
title(['spatial mean diff p=' num2str(sum(tt<t)/n)])

suptitle(sprintf('histogram of n=%d runs of random clusters \n[nonHD - HD] of corr(before,during)',n)) 

