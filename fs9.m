%room angles
% f4s2
function fs9(cellsn)
dbstop if error
load('.\data\roomdose','room');

s = {'before';'midall';'after'};tss = {'pre','dur','post'};
vs='rayleigh_score';va='rayleigh_angle';

%room figure
figure(2900);rst=0.4; st='';
clf; pr = [200 10 800 800]; set(gcf,'position',pr); 
for x=1:3 %session
    t = [cellsn.(s{x})]; ts = [t.(vs)]'; tl=arrayfun(@(z) len(z.st),t)';
    rmn=unique(room);rang= arrayfun(@(x) [t(ts>rst & room==x & tl>100).(va)],rmn,'uni',0);
    for rmi=1:len(rmn) %z to d
        subplot(3,len(rmn),(x-1)*len(rmn)+rmi);cla; 
        a=(rang{rmi});histogram(rad2deg(a),-180:20:180);xlim([-180 180]); 
        xlabel('r-angle'); axis square; [an,sd]=meanang(a);
        title(sprintf('%s rm%d n=%d a=%.0f%c%c%.0f%c',tss{x},rmn(rmi),len(a),an,176,177,sd,176));
    end
end
suptitle([st ' Rayleigh Angle by Room Number, Rayleigh Score > ' n2(rst,1)]);

end