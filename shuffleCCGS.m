
%shuffling spatial cross correlation then gridscore
function pcor0 = shuffleCCGS (c1, c2, p)
%[t1,t2] = createMsSpikeTrain(c1.st,100,c2.st); %WRONG because of aligning with xy in create ratemap
t1=createMsSpikeTrain(c1.st); nb = p.nb;
%offsets to deal with discrepency between 0 and first px py
%so that during shifting we don't put spikes before px
t2off=createMsSpikeTrain(c2.st-c2.pt(1)); %normalizes spike time relative to pt(1) = 1.
t1smooth = t1;
%t1smooth = movmean(t1,p.movmean);
rm1 = createSmoothRateMap(c1,nb,t1smooth);
%vars to speed up preallocation
mx=max(c2.px); my=max(c2.py);
rmt2 = histcounts2(c2.px,c2.py,0:mx/nb:mx,0:my/nb:my)';%save time
t2smoff = zeros(1,ceil(c2.pt(1)*1000)+len(t2off)); %zeros array, length trial offset plus time of last spike.
st2 = (1:len(t2smoff))/1000;
%shiftedi = zeros(1:len(t2off));
%size(rm1)
%delta to shift train by
d = floor(len(t2off)/p.n);%right?  %amount to shift by each step
pcor0 = zeros(1,len(p.n));
for i = 1:p.n %n+1 because wi = (i-1)*d; to include 0 shift
    %i
    %sprintf('done %.2f\%',round(i/(p.n+1),2));
    wi = (i-1)*d; %first one calculated is 0
    shiftedi = [wi+1:len(t2off) 1:wi]; %cyclic
    t2shifted = t2off(shiftedi); %shift offset spike train around by len/n
    %t2smooth = movmean(t2shifted,p.movmean);
    t2smooth = t2shifted;
    %t2smooth = movmean(t2off([wi+1:len(t2off) 1:wi]),p.movmean);
    %t2smoff = [zeros(1,ceil(c2.pt(1)*1000)) t2smooth];
    t2smoff(ceil(c2.pt(1)*1000)+1:end) = t2smooth; %all zeros from 0 until pt(1) this is ok,
    rm2 = createSmoothRateMap(c2,nb,t2smoff,rmt2,st2);%so times are always in range
    %scc = xcorr2(rm2-mean(rm2(:)),rm1-mean(rm1(:)));
    %nscc = normxcorr2(rm2,rm1); %reverse order for perspective should have negs
    scc=xcorr2(rm2,rm1);
    sig = 0.5;
    scc0 = gridscore2(xcorr2(scc),sig);
    scc1 = gridscore2(scc,sig);
    pcor0(i) = scc0;
    
    if i == 1 && scc0 == 0
       return 
    end
    
    %scc0
    %if mod(i-1,10) == 0
    
    if scc0 > 10 || i==-1
        if i==1
            figure(1);
        else
            figure(2);
        end
        colormap jet;
        subplot(221);imagesc(rm1); title(i) %axis(gca,'square');
        subplot(222);imagesc(rm2); %axis(gca,'square');
        subplot(223);imagesc(imgaussfilt(xcorr2(scc),sig)); title(scc0);
        subplot(224);imagesc(imgaussfilt(scc,sig)); title(scc1);
    end
    %end
end
end






%{

    for i = 1:length(pairs)  %   i=0;
        i = 1+i;
        c1 = cells{pairs(i,1)}.before;c2 = cells{pairs(i,2)}.before;
        t1=createMsSpikeTrain(c1.st); nb = 50;
        t1smooth = movmean(t1,p.movmean);
        rm1 = createSmoothRateMap(c1,nb,t1smooth);
        t2=createMsSpikeTrain(c2.st); nb = 50;
        t2smooth = movmean(t2,p.movmean);
        rm2 = createSmoothRateMap(c2,nb,t2smooth);
        scc = normxcorr2(rm2,rm1); %reverse order for perspective should have negs
        xcc=xcorr2(rm2,rm1);
        %imagesc(scc)
        sig = 0.5;
        s1 = gridscore(scc,sig);
        s2 = gridscore(xcc,sig);
        s3 = gridscore2(xcorr2(scc),sig);
        s4 = gridscore2(xcorr2(xcc),sig);
        tit=sprintf('%3dx%3d:%2.2f %2.2f %2.2f %2.2f \n',cells{pairs(i,1)}.ind,cells{pairs(i,2)}.ind,...,
            s1,s2,s3,s4);
        %imagesc(scc);title(tit);
        %tit
    end
%}