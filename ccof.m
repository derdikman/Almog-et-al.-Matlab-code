function [r,p,s] = ccof(a,b,rd,off)
s='';
a=a(:); b=b(:);
if nargin==4
    off=zeros(off-len(a),1); a=[a; off];b=[b; off];
    %t= a~=0 | b~=0; a=a(t);b=b(t);
end
[r, p] = corrcoef(a,b,'rows','complete');
if isnan(r(2))
    disp('nan ccof');
end
r = r(2); p = p(2);  
if nargin~=3
    rd=2;
end
f=['%0.' n2(rd) 'f'];
%s=sprintf(['r=' f ' p=' f],rnd(r,rd),rnd(p,rd));
s=sprintf(['r=' f ' %s'],rnd(r,rd),pstr(p,rd));

end

