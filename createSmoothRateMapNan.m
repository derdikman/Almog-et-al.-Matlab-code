%trs has been smoothed  - trs = movmean(createMsSpikeTrain(c.st),25)
%some vars to speed up preallocation not necessary
function rm = createSmoothRateMapNan(c,nbins,trs,rmt,st)
    %tr = createMsSpikeTrain(c.st, 100);
    %trs = movmean(tr,25);
    rm=nan(nbins);
    if len(c.st)>1
    mx = max(c.px);my = max(c.py);
    if ~exist('rmt','var')
        rmt = histcounts2(c.px,c.py,0:mx/nbins:mx,0:my/nbins:my)';
    end
    %trs is in ms. st makes each index 1ms in units of second.  
    if ~exist('st','var'); 
        %st = (1:len(trs))/1000;
        st = 1/1000:1/1000:len(trs)/1000;
    end
    
    st = st(trs>0); %to use only parts of smoothed spike train greater than 0
    
    sti = discretize(st, [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);  
    sxi = discretize(c.px(sti), 0:mx/nbins:mx);
    syi = discretize(c.py(sti), 0:my/nbins:my);
    rmss = accumarray([syi sxi], trs(trs>0), [nbins nbins]); 
    rm = rmss./rmt;
    end
        
    %rm(isnan(rm)) = 0;
    %  imgsc(rm)
end