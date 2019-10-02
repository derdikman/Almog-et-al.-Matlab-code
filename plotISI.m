function ax = plotISI(ax,c,str,p)
        bins = 10.^[-4:0.05:2];
        plot(ax,histcounts(diff(c.st),bins));%imagesc(ax, rm); %rm
        ax.XTick = round(linspace(1,length(bins),10));
        %ax.XTickLabel= round(bins(ax.XTick),4); %ADD BACK IN
        %ax.XTickLabelRotation = 270;
        %title(ax,sprintf('%s','ISI (s)'));
end

