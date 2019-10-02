function mainFrequencies()
warning off;
gsThreshB = 0.5; gsThreshM = 0; rsThreshB = 0.7;

t = hd; rw = 2; cl = 5;
for i = 1:5;%length(cells)
    c = cells{t(i)}; p = {'linewidth',0.3,'markersize',3};
    
    h= figure('position',[10 10 1600 800]); axi = [];
    %bef
    s = c.before; o=0; tr = makeMsSpikeTrain(s.st,100);
    %
    ax = subplot(rw,cl,o+1);plot(ax,s.sx,s.sy,'.');axi(end+1)=ax; %plot(ax,c.px,c.py,p);hold(ax,'on');
    ax = subplot(rw,cl,o+2);plotfft(ax,tr,1e3,500 ,5    ,1e0);axi(end+1)=ax;
    ax = subplot(rw,cl,o+3);plotfft(ax,tr,1e3,50  ,5    ,1e0);axi(end+1)=ax;
    ax = subplot(rw,cl,o+4);plotfft(ax,tr,1e3,1   ,1e-2 ,1e0-1); axis(ax,'square');
    ax = subplot(rw,cl,o+5);plotISI(ax,s,{},{}); axis(ax,'tight');
    %mid
    s = c.midall; o=cl;tr = makeMsSpikeTrain(s.st,100);
    %
    ax = subplot(rw,cl,o+1);plot(ax,s.sx,s.sy,'.');axi(end+1)=ax; %plot(ax,c.px,c.py,p);hold(ax,'on');
    ax = subplot(rw,cl,o+2);plotfft(ax,tr,1e3,500 ,5    ,1e0);axi(end+1)=ax;
    ax = subplot(rw,cl,o+3);plotfft(ax,tr,1e3,50  ,5    ,1e0);axi(end+1)=ax;
    ax = subplot(rw,cl,o+4);plotfft(ax,tr,1e3,1   ,1e-2 ,1); axis(ax,'square');
    ax = subplot(rw,cl,o+5);plotISI(ax,s,{},{}); axis(ax,'tight');
    %axi = findobj(h,'type','axes')
    for i = 1:length(axi)
        axis(axi(i),'tight');
        axis(axi(i),'square');
    end
end

 h= figure('position',[0 0 0.9 0.4],'Units','normalized');
 t = gdec; axi = []; rw = 4; cl = len(t); 
for i = 1:len(t);%length(cells)
    c = cells{t(i)}; o = rw*(i-1); p = {'parent',h,'position'};
    %bef
    sb = c.before;trb = makeMsSpikeTrain(sb.st,100);
    sm = c.midall;trm = makeMsSpikeTrain(sm.st,100);
    %
    ax = axes(p{:},apos(cl,cl*rw,i+0*cl));imagesc(ax,imgaussfilt(sb.rm,1));axi(end+1)=ax; axis(ax,'square');
    ax = axes(p{:},apos(cl,cl*rw,i+1*cl));imagesc(ax,imgaussfilt(sm.rm,1));axi(end+1)=ax; axis(ax,'square');
    
    %ax = axes(p{:},apos(cl,cl*rw,i+2*cl));plotfft(ax,trb,1e3,500 ,5,1e0);axi(end+1)=ax;
    %ax = axes(p{:},apos(cl,cl*rw,i+3*cl));plotfft(ax,trm,1e3,500 ,5,1e0);axi(end+1)=ax;
    
    ax = axes(p{:},apos(cl,cl*rw,i+2*cl));plotISI(ax,sb,{},{});axi(end+1)=ax;
    ax = axes(p{:},apos(cl,cl*rw,i+3*cl));plotISI(ax,sm,{},{});axi(end+1)=ax;
    %mid    
    %
end
    for i = 1:length(axi)
        colormap(axi(i), 'jet');
        set(axi(i),'YDir','normal');
        axis(axi(i),'tight');
        axis(axi(i),'off');
        %axis(axi(i),'square');
    end
% a.ActivePositionProperty = 'position';
% set( p, 'Widths', 600, 'Heights', 600, 'HorizontalOffsets', 100, 'VerticalOffsets', 100 )



%grid decreasing cells
t = [];
for i = 1:length(cells)
    c = cells{i};
    if c.before.gridscore > gsThreshB &&...
       c.midall.gridscore < gsThreshM &&...
       len(c.before.st) > 700 && len(c.midall.st) > 700    
        t = [t i];
    end
end; gdec = t;
%hd cells
t = [];
for i = 1:length(cells)
    c = cells{i};
    if c.before.rayleigh_score > rsThreshB
        t = [t i];
    end
end; hd = t;

%group index
cgix = [];
for i = 1:length(groups)
    g = groups{i};
    for j = 1:length(g)
        c = g(j);
        cgix(c.ind,:) = [i j] ;
    end
end
%%%
binsize = 45;
%fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_c_midall_gridscore.mat',binsize); first set
fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45);
load(fn);
[groups ~] = findSimultaneouslyRecordedCells(cells);
end