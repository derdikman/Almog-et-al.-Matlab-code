function [data,s] = apear(c1, c2)
    [s1, s2] = createMsSpikeTrain(c1.st,100,c2.st);
    p.movmean = 30;
    tcc = timeCorrelationSmoothed(s1,s2,p);
    scc = xcorr2(c2.rm-mean(c2.rm(:)),c1.rm-mean(c1.rm(:))); %reverse order for perspective
    cenr = round(size(scc,2)/2); cenc = round(size(scc,1)/2);
    scc0 = scc(cenr,cenc);
    mgs = max(c1.gridscore, c2.gridscore);
    gcc = gridscore2(xcorr2(c2.rm,c1.rm),1.4); %sigma 1.4 is good
    c1.rayleigh_score;
    data = [tcc,scc0,mgs,gcc];
    colstr = {'time cc smoothed at 30ms';
          'spatial cc at 0,0, no smoothing';
          'maximum gridscore';
          'gridscore of spatial cc';
         };

end
%{ 


update me:




%}