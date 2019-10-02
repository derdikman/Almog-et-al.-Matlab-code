 %run experiment to check corrcoef over randomized cells.
%pick cell, pick second cell from another group. 
% calculate their corr before and dur
[groups ~] = findSimultaneouslyRecordedCells(cells);
gridThreshBef = 0.3; gridThreshMid = 0.25;
[bgids, bcids] = groupsByGridThresh(groups, gridThreshBef, gridThreshMid);
gids = []; cids = [];
[gids, cids] = removeOverlappingCells(groups, bgids, bcids, -1, 0.07,[133,136]); %why 133 bad?

[s1,s2]=removeOverlappingSpikes(c1.before.st,c2.before.st,0.001);sum(s1==0) %check overlap

ccofspace= [ccof(ctsbma(:,4),ctsbma(:,5)), ccof(zb,zm)]
avspace= [mean(abs(ctsbma(:,4))), mean(abs(ctsbma(:,5))),mean(abs(zb)),mean(abs(zm)) ]
%CHECK WITH CELLS REMOVED

tb = [];  tm = [];
sb = []; sm = [];

zb = []; zm = []; p.nb = 1000;
for g = 1:len(gids)
    fprintf(' %d',g)
    for c = 1:len(cids{gids(g)})
       c1i = groups{gids(g) }(cids{gids(g )}(c )).ind; c1 = cells{c1i};
       c1i = groups{gids(g) }(cids{gids(g )}(c )).ind; c1 = cells{c1i};
       for gg= 1:len(gids)
           if gg~=g
           rc = randi(len(cids{gids(rg)}));
           rg = g;rc = c;while rc == c;rc = randi(len(cids{gids(rg)}));end; fprintf('%s','g')%SAME GROUP    
           c2i = groups{gids(rg)}(cids{gids(rg)}(rc)).ind; c2 = cells{c2i}; 
           p.movmean = 25;
           zb(end + 1) = spaceCorrelationSmoothed(c1.before,c2.before,p);
           zm(end + 1) = spaceCorrelationSmoothed(c1.midall,c2.midall,p);
           end
       end   
   end
end

%c1 = cells{22};c2 = cells{24}; c3=cells{78}; %22 b=500 in vs out

p.n=100;p.nb=500;p.movmean=25; b = []; h=[1];
for i = 1:30
    tic
g1 = randi(len(gids)); c1 = randi(len(cids{gids(g1)})); c2=c1; 
while c2==c1 || any(h==g1*100+c1*10+c2);g1 = randi(len(gids)); 
    c1 = randi(len(cids{gids(g1)}));c2 = randi(len(cids{gids(g1)}));
end;
g2 = g1; while g2==g1; g2 = randi(len(gids)); end; i
c3 = randi(len(cids{gids(g2)})); h(end+1)=g1*100+c1*10+c2;h(end+1)=g1*100+c2*10+c1; %z = [7 8 22]
c1 = groups{gids(g1)}(cids{gids(g1)}(c1));%c1 = cells{cels(z(1))}; 
c2 = groups{gids(g1)}(cids{gids(g1)}(c2));%c2 = cells{cels(z(2))};
c3 = groups{gids(g2)}(cids{gids(g2)}(c3));%c3 = cells{cels(z(3))};
%[c1.ind c2.ind c3.ind]% p
a=shuffleSpace2Correlations(c1.before,c2.before,p);b(i,1)= a(1)==max(a)|a(1)==min(a); b(i,3)=a(1);
a=shuffleSpace2Correlations(c1.midall,c2.midall,p);b(i,2)= a(1)==max(a)|a(1)==min(a); b(i,4)=a(1);
a=shuffleSpace2Correlations(c1.midall,c3.midall,p);b(i,6)= a(1)==max(a)|a(1)==min(a); b(i,5)=a(1);
toc
end; b(isnan(b(:,5)),5)=0;
stats(b(:,3)), stats(b(:,4)), stats(b(:,5))
sum(b(:,[1 2 6])) 
[ccof(b(:,3),b(:,4)) ccof(b(:,3),b(:,5))]


figure(5);  nb = 50; mm = p.movmean; colormap jet;
rm1 = createSmoothRateMap(...
    c1.before,nb,movmean(createMsSpikeTrain(c1.before.st),mm));
rm2 = createSmoothRateMap(...
    c2.before,nb,movmean(createMsSpikeTrain(c2.before.st),mm));
rm3 = createSmoothRateMap(...
    c3.before,nb,movmean(createMsSpikeTrain(c3.before.st),mm));
subplot(231);imagesc(rm1);title(sprintf('rm1 i %d %d %d',c1.ind,c2.ind,c3.ind)); axis square;
subplot(232);imagesc(rm2);title(sprintf('rm2(in) %.3f',ccof(rm2,rm1))); axis square;
subplot(233);imagesc(rm3);title(sprintf('rm3(out) %.3f',ccof(rm3,rm1))); axis square;
rm4 = createSmoothRateMap(...
    c1.midall,nb,movmean(createMsSpikeTrain(c1.midall.st),mm));
rm5 = createSmoothRateMap(...
    c2.midall,nb,movmean(createMsSpikeTrain(c2.midall.st),mm));
rm6 = createSmoothRateMap(...
    c3.midall,nb,movmean(createMsSpikeTrain(c3.midall.st),mm));
subplot(234);imagesc(rm4);title(sprintf('rm1(mid) %s','')); axis square;
subplot(235);imagesc(rm5);title(sprintf('rm2(mid) %.3f',ccof(rm5,rm4))); axis square;
subplot(236);imagesc(rm6);title(sprintf('rm3(mid) %.3f',ccof(rm6,rm4))); axis square;

%gridscore2(xcorr(rm1),1)

%OUT OF GROUP
p=[]; p.nb = 50; p.movmean=25; csb=[];csm=[];ctb=[];ctm=[];cii=[];
for i = 1:len(cels)
    fprintf(' %d',i)
    tic
    c1 = cellsn(cels(i));
    for j = i+1:len(cels)
        c2 = cellsn(cels(j));
        if ~isequal(sprintf('%s%s',c1.date,c1.id),sprintf('%s%s',c2.date,c2.id));
            csb(end+1) = spaceCorrelationSmoothed(c1.before,c2.before,p);
            csm(end+1) = spaceCorrelationSmoothed(c1.midall,c2.midall,p);
%             [t1, t2] = createMsSpikeTrain(c1.before.st, 100, c2.before.st);
%             ctb(end + 1) = timeCorrelationSmoothed(t1,t2,p);
%             [t1, t2] = createMsSpikeTrain(c1.midall.st, 100, c2.midall.st);
%             ctm(end + 1) = timeCorrelationSmoothed(t1,t2,p);
            cii(end+1,:) = [c1.ind,c2.ind];
        end
    end
    toc;
end
[mean(abs(csb)), mean(abs(csm))] %OUT OF GROUP
ccof(csb,csm)



%IN GROUP
p.nb = 100; p.movmean=25; icsb=[];icsm=[];ictb=[];ictm=[];cii=[];
for i = 1:len(cels)
    fprintf(' %d',i)
    %tic
    c1 = cells{cels(i)};
    for j = i+1:len(cels)
        c2 = cells{cels(j)};
        if isequal(sprintf('%s%s',c1.date,c1.id),sprintf('%s%s',c2.date,c2.id));
            icsb(end+1) = spaceCorrelationSmoothed(c1.before,c2.before,p);
            icsm(end+1) = spaceCorrelationSmoothed(c1.midall,c2.midall,p);
%             [t1, t2] = createMsSpikeTrain(c1.before.st, 100, c2.before.st);
%             ictb(end + 1) = timeCorrelationSmoothed(t1,t2,p);
%             [t1, t2] = createMsSpikeTrain(c1.midall.st, 100, c2.midall.st);
%             ictm(end + 1) = timeCorrelationSmoothed(t1,t2,p);
%             cii(end+1,:) = [c1.ind,c2.ind];
        end
    end
    %toc;
end
[mean(icsb) mean(icsm)] %IN GROUP



for g = 1:len(gids)
    fprintf(' %d',g)
   for c = 1:len(cids{gids(g)})
       for z= 1:100
           rg = g;while rg==g;rg = randi(len(gids));end
           rc = randi(len(cids{gids(rg)}));
           rg = g;rc = c;while rc == c;rc = randi(len(cids{gids(rg)}));end; fprintf('%s','g')%SAME GROUP          
           c1i = groups{gids(g) }(cids{gids(g )}(c )).ind; c1 = cells{c1i};
           c2i = groups{gids(rg)}(cids{gids(rg)}(rc)).ind; c2 = cells{c2i}; 
           p.movmean = 25;
           zb(end + 1) = spaceCorrelationSmoothed(c1.before,c2.before,p);
           zm(end + 1) = spaceCorrelationSmoothed(c1.midall,c2.midall,p);
       end   
   end
end
ccofspacez= [ccof(ctsbma(:,4),ctsbma(:,5)), ccof(zb,zm)]
avspacez= [mean(abs(ctsbma(:,4))), mean(abs(ctsbma(:,5))),mean(abs(zb)),mean(abs(zm)) ]
p.nb 


           %time
%            [t1, t2] = createMsSpikeTrain(c1.before.st, 100, c2.before.st);
%            tb(end + 1) = timeCorrelationSmoothed(t1,t2,p);
%            [t1, t2] = createMsSpikeTrain(c1.midall.st, 100, c2.midall.st);
%            tm(end + 1) = timeCorrelationSmoothed(t1,t2,p);
           %space
%            sb(end + 1) = spaceCorrelationSmoothed(c1.before,c2.before,p);
%            sm(end + 1) = spaceCorrelationSmoothed(c1.midall,c2.midall,p);



x = norm(eig(A)-eig(B));
norm(A-B);


a = [];
for g = 1:len(gids)
   for c = 1:len(cids{gids(g)})
       a(end+1) = groups{gids(g)}(cids{gids(g)}(c)).ind;
   end
end

a = 0;
for c = 1:len(cids)
a = a + len(cids{c});
end