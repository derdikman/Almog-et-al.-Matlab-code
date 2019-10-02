%mnxt: min and max time to make trajectory out of, set to [-1, inf]
%make empty session
% nbins = 50;
% r.ind = data.ind;
% ep.x1 = 0; ep.x2 = 0; ep.y1 = 0; ep.y2 = 0; ep.t = 1;
% es.x = 0; es.x2 = 0;  es.y = 0;  es.y2 = 0; es.ts =1;
function s = createSession(p,s,n,mnxt)
    s = createTrajectory(p.x1, p.x2, p.y1, p.y2, p.t, s.ts, mnxt); %s.x, s.x2, s.y, s.y2,
    [s.rm, s.max_r] = createRateMap(s.px, s.py, s.pt, s.sx, s.sy, s.st, true, n);
    %s.acOrig = Cross_Correlation(s.rm, s.rm); 
    s.ac = xcorr2(s.rm); 
    s.gridscore = gridscore2(s.ac, 2);
    s.module = Find_Module(imgaussfilt(s.ac, 3,'FilterDomain','spatial'));
    s.exists = false;
    if s.max_r ~= 0 %&& s.gridscore ~= -2 %&& s.max_r ~= 50 
        s.exists = true;
    end
end

