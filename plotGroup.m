function plotGroup(group, parent)
    if isempty(parent)
        parent = figure('Position', [0, 0, 2000, 1000]);
        set(gca,'LooseInset', get(gca,'TightInset'));
    end
    plotGroupPrivate(parent, group);
end

function p = fpos( ai, ncol, naxes)
    nRow = ceil( naxes / ncol ) ;
    rowH = 0.85 / nRow ;  colW = 0.85 / ncol ; %width of plot %WAS 0.75
           %offset
    colX = 0.02 + linspace( 0, 0.9, ncol+1 ) ; colX = colX(1:end-1) ; %WAS 0.05
    rowY = 0.02 + linspace( 0.9, 0, nRow+1 ) ; rowY = rowY(2:end); %0.05
    rowId = ceil( ai / ncol );
    colId = ai - (rowId - 1) * ncol;
    x = colX(colId); y = rowY(rowId);
    p = [x, y, colW, rowH];
end

function [ax ai] = axi(parent, ai, ncol, naxes)
    ax = axes('Parent',parent, 'Position', fpos( ai, ncol, naxes)) ;
    ai = ai + 1; 
end

function ai = plotTR(p);
    [ax ai] = axi(p.parent, p.ai, p.ncol, p.naxes);
    plot(ax, p.c.px, p.c.py,'k','linewidth',1.5);
    hold(ax, 'on');
    plot(ax, p.c.sx, p.c.sy, 'r.','markersize',6),...
        xlim(ax, [0 100]), ylim(ax, [0 100]);
    title(ax, p.title,'fontweight','bold');
    if ai <= p.ncol; 
        xlabel(ax, p.xl); ylabel(ax, p.yl); 
    end
end
function ai = plotRM(p)
    [ax ai] = axi(p.parent, p.ai, p.ncol, p.naxes);
    hold(ax, 'on');
    imagesc(ax, imgaussfilt(p.c.rm,2,'FilterDomain','spatial')); 
    title(ax, p.title,'fontweight','bold');
    if ai <= p.ncol; xlabel(ax, p.xl); ylabel(ax, p.yl); end
end
function ai = plotAC(p)%m,n,l
    [ax ai] = axi(p.parent, p.ai, p.ncol, p.naxes);
    imagesc(ax, imgaussfilt(p.c.ac,2,'FilterDomain','spatial'));  hold(ax, 'on');%xlim(ax, xlim); ylim(ax, ylim);%axis('xy'); axis ij; axis equal; axis off;
    x=xlim(ax); y=ylim(ax);
    %title(ax, sprintf('%sspk%d',s, length(p.sx)));
    lw = 0.5;
    %center point
    if ~isempty(p.c.module)
    if length(p.c.module.hex_peaks) == 7 %CHANGE TO EXISTS
            plot(ax, p.c.module.x, p.c.module.y,'w','LineWidth',lw);  %ADD BACK IN 
            plot(ax, p.c.module.hex_peaks(:,1),p.c.module.hex_peaks(:,2),'wo','LineWidth',lw); %ADD BACK IN
    end
    end
    title(ax, p.title,'fontweight','bold');
    xlim(ax, x); ylim(ax, y);
end
function ai = plotHD(p)
    [ax ai] = axi(p.parent, p.ai, p.ncol, p.naxes);hold(ax,'on');
    rs = [];
    c = p.r.before;rd = ft();plot(ax, 3:3:360,rd,'r');  
    c = p.r.midall;rd = ft();plot(ax, 3:3:360,rd,'b'); 
    c = p.r.after; rd = ft();plot(ax, 3:3:360,rd,'g');  
    legend(ax,{'PR','DU','PO'},'location','northeast');legend(ax,'boxoff');legend(ax,'off');
    %title(ax,sprintf('pr-r du-b po-g %.1f %.1f %.1f', rs ),'fontsize',7);
    
    %title(['\fontsize{16}black {\color{magenta}magenta ' , '\color[rgb]{0 .5 .5}teal \color{red}red} black again'],'interpreter','tex')
    function rd = ft
        rd = 0;
        if ~isnan(c.rayleigh_score)
        rs(end+1) = round(c.rayleigh_score,1);bins= 0:3:360; 
        c.si = discretize(c.st, [-Inf; mean([c.pt(2:end) c.pt(1:end-1)],2); +Inf]);
        rd = histcounts(rad2deg(c.hd(c.si))+180,bins)./(histcounts(rad2deg(c.hd)+180,bins)+0.001);
        rd = smooth(rd,15);
        end
    end
    if ai <= (p.ncol+1); xlabel(ax, p.xl); ylabel(ax, p.yl); 
    title(ax,['\fontsize{8}{\color{red}pre ','\color{blue}dur ','\color{green}pos }'],'interpreter','tex');
    end
end
function titl = plotGroupPrivate(parent, g)
    rmt = 1.5; % Rate Map Threshold Hz
    debug = '';
    ncol = 3 * 3 + 1; %3 * (1+length(g(1).middle)+1) + 1;
    naxes = length(g) * ncol;

    %MAIN LOOP
    ai = 1; p.ai = 1; p.ncol = ncol; p.naxes = naxes; p.parent = parent;
    for i = 1:length(g)
        r = g(i);
        mcross = {};
        gc = '';
        %PLOT XY
        %before
        p.c = r.before; p.title = sprintf('%sc%d:tet%d: Pre',gc, r.ind,r.tet);
        p.xl = 'position(cm)'; p.yl = 'position(cm)'; p.ai = plotTR(p);
        %middle
        if  goodCell(r.midall, rmt)
            p.c = r.midall; p.title = sprintf('%s','During');p.xl='';p.yl='';p.ai = plotTR(p);
        else;p.ai = p.ai+1;end
        %after
        if r.after.exists
            p.c = r.after; p.title = sprintf('%s','Post');p.xl='';p.yl='';p.ai = plotTR(p);
        else;p.ai = p.ai+1;end
        
        %PLOT RM
        %before
        p.c = r.before; p.title = sprintf('max bin %0.fHz',round(max(r.before.max_r(:))));
        p.xl = 'position(cm)'; p.yl = 'position(cm)'; p.ai = plotRM(p);
        %middle
        if  goodCell(r.midall, rmt)
        p.c = r.midall; p.title = sprintf('max bin %0.fHz',round(max(r.midall.max_r(:))));
        p.xl = ''; p.yl = ''; p.ai = plotRM(p);
        else;p.ai = p.ai+1;end
        %after
        if r.after.exists
        p.c = r.after; p.title = sprintf('max bin %0.fHz',round(max(r.after.max_r(:))));
        p.xl = ''; p.yl = '';p.ai = plotRM(p);
        else;p.ai = p.ai+1;end
        
        %PLOT AC
        %before
        p.c = r.before; p.title = sprintf('grid score %.1f',round(p.c.gridscore,1));
        p.ai = plotAC(p);
        %middle
        if goodCell(r.midall, rmt)
        p.c = r.midall; p.title = sprintf('grid score %.1f',round(p.c.gridscore,1));
        p.ai = plotAC(p);
        else;p.ai = p.ai+1;end
       %after
        if r.after.exists
        p.c = r.after; p.title = sprintf('grid score %.1f',round(p.c.gridscore,1));
        p.ai = plotAC(p);
        else;p.ai = p.ai+1;end
        
        %PLOT RAYLEIGH
        p.r = r;
        p.xl = 'direction'; p.yl = 'firing rate';
        p.ai = plotHD(p);
        
        %PLOT ALL CROSSED
        %{
        %[ax ai] = axi(parent, ai, ncol, naxes);%mc = allCorr(mcross,mcross); %ADD BACK 
        %mc =  padarray(mc,midmax+2-size(mc),-1,'post');mc = flip(mc); 
        %imagesc(ax, mc);caxis(ax, [-1 1]);t = mc(mc<0.99999999); %floating point madness
        %title(ax, sprintf('mn%.1f mx%.1f', round(min(t(:)),1), round(max(t(:)),1) ));
        %time correlation t2 = g(j).before.st;[p, c] = time_correlation(t1,t2);
        %}
        
        %MAKE AXIS NICE
        ax = findobj(parent,'Type','Axes');
        for i=1:length(ax)
            set(ax(i),'LooseInset', get(ax(i),'TightInset'));
            set(ax(i),'FontSize',8);
            %axis(ax(i),'off')
            set(ax(i), 'XTick', []); set(ax(i), 'YTick', []);
            colormap(ax(i),'jet')
            set(ax(i),'ydir','normal');
            axis(ax(i),'tight');
            axis(ax(i),'square')
            %axis(ax(i),'equal')
        end
        
    end
    %ADD TITLE
    %ax = axes('Parent',parent, 'Position', [.05, .95, 1, .02]); 
    ax = axes('Parent',parent, 'Position', [0.5, .98, 0, 0]);axis(ax,'off');
    txt = 'Simultaneously Recorded Entorhinal Cells [Pre, During, Post] Muscimol'; %   Trajectory   Rate Maps   Autocorrelations   Rayleigh
    text(ax, 0,0, txt,'FontSize',12, 'fontweight','bold', 'horizontalalignment','center');
    pt = {'FontSize',11, 'fontweight','bold', 'horizontalalignment','left'}; t = fpos(1,ncol,naxes);
    pa = {'Parent',parent,'Visible','off',      'Position'}; y = t(2)+t(4)+0.025;
    t = fpos(1,ncol,naxes) ; ax = axes(pa{:},[t(1) y 0 0]);text(ax, 0,0, '(a) Trajectory with Spikes',pt{:});
    t = fpos(4,ncol,naxes) ; ax = axes(pa{:},[t(1) y 0 0]);text(ax, 0,0, '(b) Firing Rate Map',pt{:});
    t = fpos(7,ncol,naxes) ; ax = axes(pa{:},[t(1) y 0 0]);text(ax, 0,0, '(c) Autocorrelation of Rate Map',pt{:});
    t = fpos(10,ncol,naxes); ax = axes(pa{:},[t(1)-.03 y 0 0]);text(ax, 0,0, '(d) Spikes by Head Direction',pt{:});
    
    
end

function g = sortByMiddle(g)
    fs = fieldnames(g);
    ie = find(strcmp(fs, 'middle'));
    a = struct2cell(g);
    s = size(a);
    a = reshape(a, s(1), []);
    a = a';
    for i = 1:length(a(:,1))
        a(i,length(fs)+1) = {length(a{i,ie})}; %middle row
    end
    a = sortrows(a, -length(fs)+1);
    a = a(:,1:length(fs));
    a = reshape(a', s);
    g = cell2struct(a, fs, 1);
end

function g = sortByEllipse(g)
    fs = fieldnames(g);
    ie = find(contains(fs, 'before'));
    a = struct2cell(g);
    s = size(a);
    a = reshape(a, s(1), []);
    a = a';
    for i = 1:length(a(:,1))
        %sort by ellipse size of *before*
        a(i,length(fs)+1) = {a{i,ie}.('module').('major_ax') * a{i,ie}.('module').('minor_ax')*pi};
        if isnan(a{i,ie}.('module').('major_ax'))  || isnan(a{i,ie}.('module').('minor_ax'))
            a(i,length(fs)+1) = {10000}; %make unknown max
        end
    end
    a = sortrows(a, length(fs)+1); %-11 for backwards
    a = a(:,1:length(fs));
    a = reshape(a', s);
    g = cell2struct(a, fs, 1);
end

function b=goodCell(c, rmt) %%LOOK AT
    if c.exists
        b = true;
        return 
    end
    disp('WOW BAD CELL PLOTGROUP');
    b= c.exists && ...
        not(isnan(c.gridscore) || ...
        c.max_r < rmt );     %...
        %c.gridscore == -2 ||  ...
        %c.max_r == 50     ||   ...
        %);
end
