function ax = plotTimeCorr(ax,c1,c2,p,str)
    title(ax, str.t);xlabel(ax, str.x);ylabel(ax, str.y); hold(ax, 'on');
    %[y, x] = compareByMovingDirection(c1,c2, params);
    %X = p{1,2}'; Y = p{1,1}';
    [t1, t2] = createMsSpikeTrain(c1.st, 100, c2.st);
    %t1 = movmean(t1,p.movmean);t2 = movmean(t2,p.movmean); p.movmean = 0;
    [Y, X] = timeCorrelationSmoothed(t1,t2,p);
    plot(ax, X, Y,'w','linewidth', 3); hold(ax, 'on');
    if ~(isfield(p,'off'))
    plot(ax,  [0 0],  [-1 1],'r','linewidth', 1); %t=ylim(ax);
    plot(ax,  0, Y(ceil(len(Y)/2)),'r.','MarkerSize', 15); %t=ylim(ax);
    text(ax, 0.5,0.9, sprintf('%.3f %.3f',round(Y(round(length(Y)/2)),3)),...
        'units','normalized','Color','r', 'Fontsize',8, 'fontweight','bold');
    end
    axis(ax, 'tight');axis(ax, 'square'); box on;
    set(ax,'Color','black');colormap(ax,'jet');
    ylim(ax, [-max(abs(Y(:))) max(abs(Y(:)))]);
    if isfield(p,'ylim')
        ylim(p.ylim);%ylim([-0.06 0.06]); %REMOVE
    end
    ax.YAxis.Exponent = 0; ax.XAxis.TickLabelFormat = '%gs';
    axis(ax,'square');
    if (isfield(p,'off'))
        set(ax, 'XTick', []); set(ax, 'YTick', []);
    end
    ylabel('corr');xlabel('lag');
end

