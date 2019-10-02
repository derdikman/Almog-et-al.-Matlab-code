function [rm, max_r] = createRateMapNan(c, nbins)
if ~isempty(c.sx)
    mx = max(c.px);
    my = max(c.py); 
    rmt = histcounts2(c.py,c.px,0:my/nbins:my,0:mx/nbins:mx); 
    rms = histcounts2(c.sy,c.sx,0:my/nbins:my,0:mx/nbins:mx);
    %add to bin where where time is 0 but there are spikes 
    rmt(rmt==0&rms~=0)=1; 
    rmt=rmt*0.02;
    rm = rms ./ rmt;
    max_r = max(rm(:));
else
    rm = nan(nbins);
    max_r = 0;
end


end

%sum(rmt(:)~=rmt2(:)) %=0 %sum(rms(:)~=rms2(:)) %=0

%pxi = discretize(px, 0:mx/nbins:mx);
%pyi = discretize(py, 0:my/nbins:my);
%rmt = accumarray([pyi, pxi], 1, [nbins nbins]);
%rmt = histcounts2(px,py,0:mx/nbins:mx,0:my/nbins:my)'; %sum(rmt(:)~=rmt2(:)) 0=0

%rmt(rmt<1) = 1e6; %NO NANS??
%sxi = discretize(sx, 0:mx/nbins:mx); %px(si)
%syi = discretize(sy, 0:my/nbins:my); %py(si)
%rms = accumarray([syi,sxi], 1, [nbins nbins]);
%rms = histcounts2(sx,sy,0:mx/nbins:mx,0:my/nbins:my)'; % sum(rms(:)~=rms2(:))% =0

%t = rms-rmt; t(t<0)=0; rmt = rmt + t; 