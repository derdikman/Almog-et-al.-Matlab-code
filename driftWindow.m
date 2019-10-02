function [ windowed ] = driftWindow(c, winsecs)
s = c;
winend = s.pt(1) + winsecs;
winref.x = s.px(1); winref.y = s.py(1);    
windowed{1}.x = []; windowed{1}.y = [];
windowed{1}.start = s.pt(1); windowed{1}.end = s.pt(1);
windowed{1}.sx = []; windowed{1}.sy = []; 
windowed{1}.t = []; windowed{1}.pt = []; 
spki = 1; wini = 1;
%skips last window
minrms = inf;
for i = 1:length(s.pt)
    if s.pt(i) >= winend || i == length(s.pt)
        %collect
        windowed{wini}.sx = [0];
        windowed{wini}.sy = [0];
        windowed{wini}.st = [1];
        t = 1;
        while spki <= length(s.st) &&  s.st(spki) < winend
            windowed{wini}.sx(t) = s.sx(spki) - winref.x; 
            windowed{wini}.sy(t) = s.sy(spki) - winref.y;
            windowed{wini}.st(t) = s.st(spki);
            spki = spki + 1; t = t + 1;
        end
        windowed{wini}.end = s.pt(i);
        t = windowed{wini};
        t.x = t.x-min(t.x); t.y = t.y-min(t.y); t.sx = t.sx-min(t.sx); t.sy = t.sy-min(t.sy);
        windowed{wini}.rm = createRateMap(t.x, t.y, t.t, t.sx, t.sy, t.st,true,50);
        minrms = min(minrms, min(size(windowed{wini}.rm)));
        %reset
        winend = s.pt(i) + winsecs;
        winref.x = s.px(i); winref.y = s.py(i); %sets reference to beginning of window
        %initialize next window
        wini = wini+1;
        windowed{wini}.x = []; windowed{wini}.y = [];
        windowed{wini}.sx = []; windowed{wini}.sy = [];
        windowed{wini}.t = []; windowed{wini}.st = [];
        windowed{wini}.start = s.pt(i);
    end
    windowed{wini}.t(end+1) = s.pt(i);
    windowed{wini}.x(end+1) = s.px(i) - winref.x;
    windowed{wini}.y(end+1) = s.py(i) - winref.y;
end
windowed = windowed(1:end-1); %last window incomplete so cutoff

function plotDrift()
    r = length(windowed) + 1;
    c = 3;
    h = m.fig;
    tot.rm = zeros(minrms);
    for i = r:r
        %set(m.parent,'LooseInset', get(m.parent,'TightInset'));
        a = windowed{i-1};
        %subplot(r,c,c*(i-1)+1);
        ax = subplot(131, 'Parent', m.parent);
        plot(ax, a.x,a.y); hold(ax, 'on'); scatter(ax, a.sx,a.sy,'r.'); axis(ax,'square');
        ax = subplot(132, 'Parent', m.parent);
        %subplot(r,c,c*(i-1)+2);
        imagesc(ax, a.rm); set(ax,'ydir','normal'); axis(ax, 'square');
        myColorMap = jet(256); myColorMap(1,:) = 1; colormap(ax, myColorMap);
        %subplot(r,c,c*(i-1)+3);
        ax = subplot(133, 'Parent', m.parent);
        a.ac = Cross_Correlation(a.rm, a.rm);
        imagesc(ax, a.ac); set(ax,'ydir','normal'); axis(ax,'square');
        colormap(ax,'jet');
        titl = sprintf('win %ds start %.0fs.png',winsecs, a.start);
        tot.rm = tot.rm + a.rm(1:minrms,1:minrms);
    end
end
%{
g = figure(); colormap jet;
subplot(121);
imagesc(m.rm); set(gca,'ydir','normal'); axis square;
subplot(122);
m.ac = Cross_Correlation(m.rm, m.rm);
imagesc(m.ac); set(gca,'ydir','normal'); axis square;
titl = sprintf('mean rate map win %ds.png',winsecs);
%}
%print(gcf, 'test.png', '-dpng', '-r300' );   %save file as PNG w/ 300dpi
disp('')
end

