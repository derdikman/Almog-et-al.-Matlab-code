%trs has been smoothed  - trs = movmean(createMsSpikeTrain(c.st),25)
%some vars to speed up preallocation not necessary
function rm = createSmoothRateMap(c,nbins,trs,rmt,st)
    %tr = createMsSpikeTrain(c.st, 100);
    %trs = movmean(tr,25);
     mx = max(c.px);my = max(c.py);
%     pxi = discretize(c.px, 0:mx/nbins:mx);
%     pyi = discretize(c.py, 0:my/nbins:my);
%     rmt = accumarray([pyi pxi], 1, [nbins nbins]);toc %change to histcounts
    if ~exist('rmt','var')
        rmt = histcounts2(c.px,c.py, 0:mx/nbins:mx,0:my/nbins:my)'; 
        %rmt(:)=1;%REMOVE
        %rmt2 = histcounts2(c.px,c.py, [nbins nbins])'; irregular bins
    end
    if ~exist('st','var');st = (1:len(trs))/1000;end
    %rmt = accumarray([pyi pxi], t, [nbins nbins]);
    %rmt(rmt<1) = 1e10; Taken care of by isnan() below
    trs0 = trs>0; st0 = st(trs0); %to use only parts of smoothed spike train greater than 0
    %sti = discretize(st, [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);
    sti = discretize(st0, [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);  
    sxi = discretize(c.px(sti), 0:mx/nbins:mx);
    syi = discretize(c.py(sti), 0:my/nbins:my);
    %rmss = accumarray([syi sxi], trs, [nbins nbins]);
    rmss = accumarray([syi sxi], trs(trs0), [nbins nbins]); %USE THIS ONE
    rm = rmss./rmt;
    
    rm(isnan(rm)) = 0; 
    
    t = rm(:);td=std(t(t>0));n=0.5;    %stats(t(t>0));
    %rm(rm<n*td) = 0;
    %rm(rm>0)=1;
    %rm(rmss>0)=1;%REMOVE
    %rm(rmss<10)=0;%REMOVE
end

%MATLAB docs: "The behavior of discretize is similar 
%to that of the histcounts function. Use histcounts to find 
%the number of elements in each bin. On the other hand, use discretize to find which bin each element belongs to (without counting)."


%{
figure;imagesc(rm);
figure;imagesc(c.rm);
figure;imagesc(rmt);
figure;imagesc(rmss);
figure;plot(c.sx,c.sy,'.');
figure;histogram(sxi,50)
figure;histogram(c.px,50)
figure;histogram(st,50)
figure;histogram(sti,50)



%}