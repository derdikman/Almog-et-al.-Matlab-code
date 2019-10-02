%mnxt: min and max time to make trajectory out of
function c = createTrajectory(px1, px2, py1, py2, pt, st, mnxt) %sx1, sx2, sy1, sy2,
    pi = mnxt(1) <= pt & pt <= mnxt(2); si = mnxt(1) <= st & st<=mnxt(2);
    pt = double(toCol(pt(pi))); st = double(toCol(st(si)));
    minx = min(min(px1),min(px2));miny = min(min(py1),min(py2));
    px1 = double(toCol(px1(pi)) - minx +.00001); 
    py1 = double(toCol(py1(pi)) - minx +.00001);
    px2 = double(toCol(px2(pi)) - miny +.00001); 
    py2 = double(toCol(py2(pi)) - miny +.00001);
    
    %PATCH TRAJECTORY
    dt = median(diff(pt));
    p1 = patchTrajectoryLinear(pt,px1,py1,dt,1.1*dt); 
    p2 = patchTrajectoryLinear(pt,px2,py2,dt,1.1*dt); 
    c.pt = p1.t; %PT%
    
    %case for no spikes
    %c.st = st; REMOVED, NEED BACK IN??
    if isempty(st)
        st = c.pt(1);
    end
    
    %get closest time in px
    si = 1; 
    if length(c.pt) > 1
        si = discretize(st, [-Inf; mean([c.pt(1:end-1) c.pt(2:end)],2); +Inf]);
    end
    
    %base sx sy on closest time of spike to data in px py
    winms = 10;% smoothing window for rayleigh score 
    [c.rayleigh_score, c.rayleigh_angle, c.hd c.shd] =...
        rayleigh_score(c.pt,p1.x,p1.y,p2.x, p2.y,...
        p1.x(si),p1.y(si),p2.x(si),p2.y(si),winms);
    
    c.px =  mean([p1.x, p2.x], 2);
    c.px = c.px - min(c.px) + 0.00001; %no zeros
    c.py =  mean([p1.y, p2.y], 2);
    c.py = c.py - min(c.py) + 0.00001;
    %base sx sy on closest time of spike to data in px py
    c.st = st;
    c.sx = c.px(si);
    c.sy = c.py(si);
    
%     c.sx = mean([sx1, sx2], 2); c.sx = c.sx - min(c.sx) + 0.00001;
%     c.sy = mean([sy2, sy2], 2); c.sy = c.sy - min(c.sy) + 0.00001;
    %{figure();plot(c.px, c.py);hold on; plot(c.sx, c.sy,'.');%}
end

