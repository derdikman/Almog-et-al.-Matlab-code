function ax = plotAC(ax,c,str)
    title(ax, str.t);xlabel(ax, str.x);ylabel(ax, str.y); hold(ax, 'on');
    imagesc( imgaussfilt(c.ac,2,'FilterDomain','spatial')); lw = 0.5;
    %title(ax, sprintf('%sspk%d',s, length(p.sx)));
    %center point ADD BACK
%     if length(c.module.hex_peaks) == 7 %CHANGE TO EXISTS
%             plot(ax, c.module.x, c.module.y,'w','LineWidth',lw); 
%             plot(ax, c.module.hex_peaks(:,1),c.module.hex_peaks(:,2),'wo','LineWidth',lw);
%     end
    xlabel('cm');ylabel('cm'); ax=gca;
    set(ax,'ydir','normal');colormap(ax,'jet');axis(ax,'square');
    if exist('str','var') && (isfield(str,'lim'))
       tmax=round(str.lim/10)*10;
    else
       tmax=round(max([max(c.px) max(c.py)])/10)*10;
    end
    %xmax=round(max(c.px)/10)*10;ymax=floor(max(c.py)/10)*10;
    t=slim(ax);t=floor(t/10)*10;xlim(t);ylim(t);
    t(3)=t(end);t(2)=round(t(end)/2);d=4; t(1)=t(1)+d; t(end)=t(end)-d; 
    xticks(t);yticks(t); t={n2(-tmax);'0';n2(tmax)};
    xticklabels(t);yticklabels(t);box off;
end
%   figure(991)
