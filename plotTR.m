function ax = plotTR(ax,c,str,p)%SETS LIMITS OF ARENA TO 100
    title(ax, str.t);xlabel(ax, str.x);ylabel(ax, str.y);
    plot(ax, c.px, c.py,'k','linewidth',1.5); hold(ax, 'on');
    plot(ax, c.sx, c.sy, 'r.','markersize',2.7);% 4 , 2.7
    set(ax,'ydir','normal');colormap(ax,'jet');axis(ax,'square');
    xlabel('cm');ylabel('cm');ax=gca;
    if exist('str','var') && (isfield(str,'lim'))
       t=[0 str.lim];
    else
       t=[0 max([max(c.px) max(c.py)])];
    end
    t=round(t/10)*10;xlim(t);ylim(t);
    d=2; tt(1)=t(1)+d; tt(2)=t(end)-d; 
    xticks(tt);yticks(tt); tt={}; tt{1}=n2(t(1));tt{2}=n2(t(2));
    xticklabels(tt);yticklabels(tt);
    box off;
    %set(ax,p);
end

