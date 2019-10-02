%returns all values from start to end.
function cw=windowsesh(c,i1,i2, st, et)
    if exist('st','var')
        i1 = find(c.pt>=st,1,'first'); %[~,i1] = max(c.pt>=st)
        i2 = find(c.pt<=et,1,'last');%[~,i2] = max(c.pt(c.pt<=et))
    end
    if isempty(i1) || isempty(i2) 
        cw.pt = 1; cw.px=0;cw.py=0;
        cw.st = 1; cw.sx=0;cw.sy=0;
        cw=[];%?? why have previous lines?
    else
        cw.pt=c.pt(i1:i2);
        cw.px=c.px(i1:i2);
        cw.py=c.py(i1:i2);
        cw.hd=c.hd(i1:i2);
        cw.st=c.st(c.st>=cw.pt(1) & c.st<=cw.pt(end));
        cw.sx=c.sx(c.st>=cw.pt(1) & c.st<=cw.pt(end));
        cw.sy=c.sy(c.st>=cw.pt(1) & c.st<=cw.pt(end));
    end
end


