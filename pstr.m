function s=pstr(p,r)
% p=0.004; r=2; 
%will be mathemetically wrong for some instances but this is what was asked for

    p=rnd(p,r)  ;
    if p == 0
        s = ['p<' n2(1/10^(r+1),r+1)];
    else
        s=['p=' n2(rnd(p,r),r)];
    end
   
end