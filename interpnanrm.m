%interpolates nans in ratemap
function rmintrp=intrepnanrm(rm)
    rmo=rm;
    [m n]=size(rm);
    [x y]=meshgrid(1:n,1:m);
    x=x(:);y=y(:);rm=rm(:);
    nni=isnan(rm); nx=x;ny=y;
    nx(nni)=[];ny(nni)=[];rm(nni)=[];
    f = scatteredInterpolant(nx,ny,rm,'natural','nearest');
    rmintrp=f(x,y);
    if len(rmintrp) ~= len(rmo(:));
        'error cannot interp ratemap'
        rmintrp = zeros(m,n); rmintrp(1,1)=1;
    else
        rmintrp=reshape(rmintrp,m,n);
    end
end


%{
rm=cellsn(2).before.rm;
%rm=rm(25:50,:)
figure(1);
subplot(211);
rmt=rm; %rmt(isnan(rmt))=99;
imgsc(rmt)
subplot(212);
imgsc(interpnanrm(rm),3);
%}
