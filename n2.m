function s=n2(n,rd)
    r=2;
    if nargin==2
        r=rd;
    end
    s=sprintf(['%0.' num2str(r) 'f' ],round(n,r));
    if floor(n)==n                      
     s = num2str(n);
    end
end