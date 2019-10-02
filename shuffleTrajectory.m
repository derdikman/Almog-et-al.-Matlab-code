function [ output_args ] = shuffleTrajectory( input_args )
    
    
    pairs = [2,7;11,12;11,13;11,14;12,13;12,14;13,14;22,24;22,28;24,28;34,35;34,37;34,39;34,41;35,37;35,39;35,41;37,39;37,41;39,41;51,53;61,62;61,66;61,68;62,66;62,68;66,68;77,78;77,79;78,79;80,81;82,86;100,101;100,103;100,104;100,105;100,111;101,103;101,104;101,105;101,111;103,104;103,105;103,111;104,105;104,111;105,111;112,114;112,115;112,116;112,117;112,118;112,119;114,115;114,116;114,117;114,118;114,119;115,116;115,117;115,118;115,119;116,117;116,118;116,119;117,118;117,119;118,119;171,173;171,174;171,176;171,177;173,174;173,176;173,177;174,176;174,177;176,177;178,179;178,180;178,182;178,184;179,180;179,182;179,184;180,182;180,184;182,184;188,189;199,201;199,209;201,209;230,231;274,276;274,278;276,278;282,283;282,287;283,287];
    
    b = []; m = []; a = [];
    
    figure(34);c = cells{2}.before;tr = createMsSpikeTrain(c.st,100);
    for i = 1:20%length(pairs{
        %c = cells{pairs(i,1)}; c = c.before;
        %tr = makeMsSpikeTrain(c.st,100);
        win = 2^(i-1);
        trs = movmean(tr,win);
        st = [1:len(trs)]/1000;
        nbins = 50;
        mx = max(c.px);
        my = max(c.py);
        sti = discretize(st, [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);
        sxi = discretize(c.px(sti), 0:mx/nbins:mx);
        syi = discretize(c.py(sti), 0:my/nbins:my);
        rmss = accumarray([syi sxi], trs, [nbins nbins]);
        pxi = discretize(c.px, 0:mx/nbins:mx);
        pyi = discretize(c.py, 0:my/nbins:my);
        %t = diff(c.pt(pw)); t = [median(t); t];
        rmt = accumarray([pyi pxi], 1, [nbins nbins]);
        %rmt = accumarray([pyi pxi], t, [nbins nbins]);
        rmt(rmt<1) = 1e10;
        rm = rmss./rmt;
        subplot(5,4,i);
        %imagesc(imgaussfilt(rm,1));colormap jet; axis(gca,'square');
        imagesc(rm);colormap jet; axis(gca,'square');
        %plot(trs)
        set(gca,'YDir','normal'); axis(gca,'tight');
        title(sprintf('win %.dms',win),'color','m','fontsize', 12);
    end 
    
end

