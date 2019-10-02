function pcor0 = spaceCorrelationSmoothed(c1,c2,p)
    %     rm1
    t1=createMsSpikeTrain(c1.st); %nb = p.nb;
    t1s = movmean(t1,p.movmean);
    %rm1 = createSmoothRateMap(c1,nb,t1s);
    rm1 = createSmoothRateMapNan(c1,p.nb,t1s);
    %     rm2
    t2=createMsSpikeTrain(c2.st);
    t2s= movmean(t2,p.movmean);
    %rm2 = createSmoothRateMap(c2,nb,t2s);
    rm2 = createSmoothRateMapNan(c2,p.nb,t2s);

    pcor0=ccof(rm2,rm1);
    
%     figure()
%     subplot(211); imagesc(rm1);
%     subplot(212); imagesc(rm2);

end