function ax = plotCCTime(c1,c2,str,ax,p)%SETS LIMITS OF ARENA TO 100
    title(ax, str.tt);xlabel(ax, str.xt);ylabel(ax, str.yt);hold(ax, 'on');
    
    [p, co] = compareByMovingDirection(c1,c2, p);
    X = p{1,2}';
    if p.nbins == 1 %Time correlation, no angle
        Y = p{1,1}';
        if false%sigma ~= 0
            gwindow = round(length(Y)/10); %gaussian win: 15 vs 30 very similar default matlab is 5;
            Y = smoothts(Y,'g',gwindow,10^sigma); %DELETE??
        end
        plot(ax, X, Y,'linewidth', 3); hold(ax, 'on');
        plot(ax,  [0 0],  [-1 1],'r','linewidth', 1);
        Y = p{1,1}';
        %plot at halfway
        text(ax, 0,0, sprintf('%.3f %.3f', round(Y(round(length(Y)/2)),3)), 'Color','r');
        axis(ax, 'tight'); 
        ax.Color = 'black';
        
    else %plot by degree bins
        Y = reshape(cell2mat(p(:,1)),[],nbins);
        %[~,ax,AX] = plotmatrix(X,Y,'-');
        Y = Y'; %Y2 = Y' %imagesc(imgaussfilt(Y2, 2));set(gca,'ydir','normal');
        if sigma ~= 0
            Ysm = imgaussfilt([Y;Y;Y], sigma,'FilterDomain','spatial'); %for nans
            Ysm = Ysm(size(Y,1)+1:size(Y,1)*2,1:size(Y,2));
            if any(~isnan(Ysm(:))) %if any are not nan
                Y = Ysm;
            end
        end
        imagesc(ax, Y); cmax = max(abs(Y(:))); colormap(ax,'jet');
        set(ax,'ydir','normal'); caxis(ax, [-cmax cmax]);
        axis(ax,'tight')
        set(ax,'xtick',[1 round(length(X)/2) length(X)]);
        set(ax,'xticklabel',[X(1) 0 X(end)]);%for imagesc
        set(ax,'ytick',1:length(p(:,3)));
        set(ax,'yticklabel', mean(cell2mat(p(:,3))')); %#ok<*UDIM>
        if c1i==1 && c2i==2 %used to be AX)
            set(ax,'fontweight','bold');
            ylabel(ax,sprintf('MD bin (b=%.1f°)',360/nbins));
            xlabel(ax,sprintf('%s sbin%.2fs','Xcorr lag(s)', ...
                params.time_bin_secs));
        end
    end
    
    axis(ax,'tight'); set(ax,'ydir','normal');colormap(ax,'jet');
end