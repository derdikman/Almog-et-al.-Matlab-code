function stat = stats( i,nans,infs)
    i = toCol(i);
    if nargin > 1 
        i=i(~isnan(i));
        if nargin==3; i=i(i~=inf&i~=-inf);end;
    end
    a.n = len(i);
    a.mean = mean(i);
    a.median=median(i);
    a.min=min(i);
    a.max=max(i);
    a.std=std(i);
    stat=a;
    stat
    %datastats(i)
end

