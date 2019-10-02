function ax = plotHD(ax,c,str,w180)%SETS LIMITS OF ARENA TO 100
    title(ax, str.t);xlabel(ax, str.x);ylabel(ax, str.y);hold(ax, 'on');
    l=[0 360];
    if nargin==4;l=[-180 180];end
    bins= l(1):1:l(2);
    c.si = discretize(c.st, [-Inf; mean([c.pt(2:end) c.pt(1:end-1)],2); +Inf]);    
    rd = histcounts(wa(c.hd(c.si),l),bins)./(histcounts(wa(c.hd,l),bins)+0.00001);
    rd = movmean(rd,round(len(rd)/5));
    p=plot(ax,bins(2:end),rd/1); p.LineWidth = 1.5; 
    xlabel('\theta'); ylabel(' ');
    set(ax,'xtick',l,'yticklabel',{},'ylim',[0 max(rd)],... 
        'xticklabel',{sprintf('%d%c',l(1),char(176)); sprintf('%d%c',l(2),char(176))})
           
    axis(ax, 'square');colormap(ax,'jet'); 
    %figure, plot(rd/sum(rd)), ylim([0 1/10000])
    
    function o=wa(a,l)
        if l(1)==0; o= wrapTo360(rad2deg(a));
        else; o=rad2deg(a); end
    end
end

