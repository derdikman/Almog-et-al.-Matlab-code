function ax = plotSpaceCorr(ax,c1,c2,p,str,cmeth)%SETS LIMITS OF ARENA TO 100
    title(ax, str.t);xlabel(ax, str.x);
    ylabel(ax, str.y);hold(ax, 'on');
    cmet='g';
    if nargin==6
        cmet=cmeth;
    end
    if isequal(cmet,'n')
        cc = xcorr2n(c2.rm,c1.rm); %reverse order for perspective
    else
        cc = xcorr2g(c2.rm,c1.rm); %reverse order for perspective
    end
    
    cc = imgaussfilt(cc,2);
    cenr = round(size(cc,2)/2); cenc = round(size(cc,1)/2);
    imagesc(ax, cc);
    %text(cenr,cenc,sprintf('%.1f',cc(cenr,cenc)),...'Color','w');%'FontSize',10);
    %plot(ax, size(cc,2)/2,size(cc,1)/2,'md','MarkerFaceColor','w','MarkerSize',7); %mark 0,0
    %intersection union
    if isfield(p,'module')
        plot(ax, c1.module.x, c1.module.y,'w');
        plot(ax, c2.module.x, c2.module.y,'r');
        title(ax, sprintf('Int/Union=%.2f',Intersection_Over_Union(...
            c2.module.x,c2.module.y,c1.module.x,c1.module.y) ),'fontsize',9);
    end

    %axis(ax,'off');axis(ax,'tight');
    axis(ax,'square'); set(ax,'ydir','normal');
    colormap(ax,'jet');
    if (isfield(p,'off'))
        set(ax, 'XTick', []); set(ax, 'YTick', []);
    end
    xlabel('cm');ylabel('cm'); ax=gca;
    if exist('str','var') && (isfield(str,'lim'))
        tmax=round(str.lim/10)*10;
    else
        tmax=round(max([max(c1.px) max(c1.py)])/10)*10;
    end
    %xmax=round(max(c.px)/10)*10;ymax=floor(max(c.py)/10)*10;
    t=slim(ax);t=floor(t/10)*10;xlim(t);ylim(t);
    t(3)=t(end);t(2)=round(t(end)/2);d=4; t(1)=t(1)+d; t(end)=t(end)-d;
    xticks(t);yticks(t); t={n2(-tmax);'0';n2(tmax)};
    xticklabels(t);yticklabels(t);box off;
end
