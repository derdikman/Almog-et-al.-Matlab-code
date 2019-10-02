function ctsbma=tsbma(cellsn,pairs)
ctsbma=[];
ss={'before' 'midall','after'};
mvmn=25; nb = 50;
i=0;
for p=pairs'
    i=i+1;
    j=0;
    for s= [ss(:)]'
      j=j+1;
      ct=nan;cs=nan;
      a=cellsn(p(1)).(s{:});b=cellsn(p(2)).(s{:});
      if len(a.st)>1 && len(b.st)>1
      [at, bt] = createMsSpikeTrain(a.st, -1, b.st);
      at = movmean(at,mvmn);arm=createSmoothRateMapNan(a,nb,at);
      bt = movmean(bt,mvmn);brm=createSmoothRateMapNan(b,nb,bt);     
      ct = ccof(at, bt); %off=sum(at==0);
      cs = ccof(arm,brm,2,len(at)); %check order
      end
      ctsbma(i,j)=ct;ctsbma(i,j+3)=cs;
    end   
end


end