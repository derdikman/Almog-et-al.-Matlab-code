function [h, vers] =  plotByDirectionMainTimeTab(params)
    %PARAMETERS
    warning off;
    params.gid
    binl = params.binspike;%0.06 0.1;%.3
    lag = params.lag*1000;
    params.lag = params.lag*1000;
    nbins = round(360/params.bindeg); %rounds good!  fprintf('nbins %d\n', nbins);
    c = 2;
    
    sigma = params.sigma
    if sigma ~= 0
        sigma = round(sigma)
        %%sigma =  (c^(params.sigma/10)-1)/(c^10-0.99999) %999 good
        %sigma =  (log(params.sigma))/(log(100))-1e-8 %999 bad
        %sigma = params.sigma/100-1e-6 %bad
    end
    params.sigma = sigma;
    v = params.sesh;
    params.number_degree_bins = nbins; params.lag_max_secs = round(lag); params.time_bin_secs = binl;
    sesh = params.sesh;
    good = params.good;
    if isa(params,'App')
        h = params.UIFigure;
        parent = params.parent;
    elseif isfield(params,'parent')
        parent = params.parent;
        h = params.fig;
    else
        h = figure('Position', [0, 0, 2000, 1000]);
        set(gca,'LooseInset', get(gca,'TightInset')); 
        parent = h; 
    end
        
    % plot each cell %
    c = length(good)+1; r = length(good);  if r == 1; r=3; end; z = 0;
    if ~isempty(good) && length(good) >= 2        
        %wt = parent.Position(3); dwt = floor(wt/c);ht = parent.Position(4); dht = floor(ht/c);
        %for each cell, for each direction x each other cell
        mxcrs = [];
        corrAxes = [];
        for c1i = 1:length(good)-1
            for c2i = c1i+1:length(good)
                z = z + 1;
                if strcmp(params.sesh,'before')
                    c1 = good(c1i).before; c2 = good(c2i).before; %clast = good(r).before;
                elseif strcmp(params.sesh,'midall')
                    c1 = good(c1i).midall; c2 = good(c2i).midall;   %clast = good(r).midall;
                elseif strcmp(params.sesh,'after')
                    c1 = good(c1i).after; c2 = good(c2i).after;   %clast = good(r).after;
                end 
                c1.ind = good(c1i).ind; c2.ind = good(c2i).ind;  %clast.ind = good(r).ind;
                t1 = length(c1.st); t2 = length(c2.st);
                %remove overlapping spikes
                rm1 = c1.rm; rm2 = c2.rm;
                [s1i, s2i] = removeOverlappingSpikes(c1.st,c2.st,params.overlap);
                c1.sx = c1.sx(s1i); c2.sx = c2.sx(s2i); c1.sy = c1.sy(s1i); c2.sy = c2.sy(s2i); c1.st = c1.st(s1i); c2.st = c2.st(s2i);
                fprintf('pre spikes removed %dx%d = %d     t1 %.1f t2 %.1f\n',c1.ind,c2.ind,t1-length(c1.st),length(c1.st)/t1,length(c2.st)/t2);
                c1.rm = createRateMap(c1.px, c1.py, c1.pt, c1.sx, c1.sy, c1.st, true, 50); 
                c2.rm = createRateMap(c2.px, c2.py, c2.pt, c2.sx, c2.sy, c2.st, true, 50);
                c1.ac = xcorr2(c1.rm); c2.ac = xcorr2(c2.rm); c1.gridscore = gridscore2(c1.ac,3); c2.gridscore = gridscore2(c2.ac,3);
                c1.module = Find_Module(imgaussfilt(c1.ac,3));c2.module = Find_Module(imgaussfilt(c2.ac,3));
                
                %COMPUTE CORRELATions
                [p co] = compareByMovingDirection(c1,c2, params);
                
                %CORRELATION PLOT
                az = subplot(r, c, c*(c1i-1) + c2i-1,'Parent', parent);
                X = p{1,2}';
                if nbins == 1 %Time correlation, no angle
                    Y = p{1,1}';
                    imax = p{1,3}; rawmax = imax(1); 
                    if false%sigma ~= 0
                       gwindow = round(length(Y)/10); %gaussian win: 15 vs 30 very similar default matlab is 5;
                        Y = smoothts(Y,'g',gwindow,10^params.sigma); %DELETE??
                    end
                    plotmax = max(abs(Y(:)));
                    plot(az, X, Y,'linewidth', 3); hold(az, 'on'); 
                    plot(az,  [0 0],  [-1 1],'r','linewidth', 1);
                    Y = p{1,1}';
                    %plot at halfway
                    text(az, 0,0, sprintf('%.3f %.3f', round(Y(round(length(Y)/2)),3)), 'Color','r',...
                        'fontsize',16); 
                    axis(az, 'tight'); %ylim(az, [-1 1]); %lim done at end
                    axis(az, 'square');
                    az.Color = 'black'; 
                    
                    
                else %plot by degree bins
                    Y = reshape(cell2mat(p(:,1)),[],nbins);
                    %[~,ax,AX] = plotmatrix(X,Y,'-');
                    Y = Y'; %Y2 = Y' %imagesc(imgaussfilt(Y2, 2));set(gca,'ydir','normal');
                    rawmax = max(abs(Y(:)));
                    plotmax = max(abs(Y(:)));
                    if params.sigma ~= 0
                        Ysm = imgaussfilt([Y;Y;Y], params.sigma,'FilterDomain','spatial'); %for nans
                        Ysm = Ysm(size(Y,1)+1:size(Y,1)*2,1:size(Y,2));
                        if any(~isnan(Ysm(:))) %if any are not nan
                            Y = Ysm;
                        end    
                    end
                    imagesc(az, Y); cmax = max(abs(Y(:))); colormap(az,'jet'); 
                    set(az,'ydir','normal'); caxis(az, [-cmax cmax]);
                    axis(az,'tight')
                    set(az,'xtick',[1 round(length(X)/2) length(X)]);
                    set(az,'xticklabel',[X(1) 0 X(end)]);%for imagesc
                    set(az,'ytick',1:length(p(:,3)));
                    set(az,'yticklabel', mean(cell2mat(p(:,3))')); %#ok<*UDIM>
                    if c1i==1 && c2i==2 %used to be AX)
                        set(az,'fontweight','bold');
                        ylabel(az,sprintf('MD bin (b=%.1f°)',360/nbins));
                        xlabel(az,sprintf('%s sbin%.2fs','Xcorr lag(s)', ...
                        params.time_bin_secs));
                    end
                end
                mxcrs = [mxcrs plotmax];
                corrAxes = [corrAxes az];
                %title(az, sprintf('%dx%d|mx|%.2f',c1.ind,c2.ind,rawmax)); %WITH MAX CORR
                title(az, sprintf('%dx%d',c1.ind,c2.ind),'fontweight','bold');
                set(az,'FontSize',9);
                %loop for plotmatrix(?) for i=1:length(ax) %title(ax(i),{'Very Nice'}) axis(ax(i),'off');     
                 %end plot correlation
                 
                %PLOT CROSS SPATIAL
                az = subplot(r, c, r*c-(c*(c1i-1) + c2i),'Parent', parent); hold(az,'on');
                cc = xcorr2(c2.rm,c1.rm); %reverse order for perspective
                imagesc(az, cc); xlabel(az, sprintf('c%d x c%d',c1.ind,c2.ind),'fontweight','bold');colormap(az,'jet');
                plot(az, size(cc,2)/2,size(cc,1)/2,'md','MarkerFaceColor','w','MarkerSize',7);
                title(az, sprintf('%dx%d',c1.ind,c2.ind),'fontweight','bold');
                %intersection union
                %plot(az, c1.module.x, c1.module.y,'w'); plot(az, c2.module.x, c2.module.y,'r');  %PUT BACK
%                xlabel(az, sprintf('Int/Union=%.2f',Intersection_Over_Union(...,MISSING TOOLBOX
 %                   c2.module.x,c2.module.y,c1.module.x,c1.module.y) ),'fontsize',9);
                axis(az,'equal'); axis(az,'tight'); set(az,'ydir','normal','xticklabel',[],'yticklabel',[]);
            end
            
            %PLOT RM
            az = subplot(r,c,c1i*c-1,'Parent', parent); imagesc(az, imgaussfilt(c1.rm,1));
            title(az, sprintf('c%d',c1.ind));
            axis(az,'equal'); axis(az,'tight'); set(az,'ydir','normal','xticklabel',[],'yticklabel',[]);
            xlabel(az,sprintf('Max %.0fHz',c1.max_r));
            colormap(az,'jet');
            
            % PLOT AC
            az = subplot(r,c,c1i*c,'Parent', parent); imagesc(az, imgaussfilt(c1.ac,2));
            hold(az,'on'); %plot(az, c1.module.x, c1.module.y,'w');
            title(az, sprintf('gsc %.1f',round(c1.gridscore,1)));
            axis(az,'equal'); axis(az,'tight'); set(az,'ydir','normal','xticklabel',[],'yticklabel',[]);
            colormap(az,'jet');
        end
        
        %ADJUST XCORR LIMS
        mxcr = max(abs(mxcrs));
        if nbins == 1 
            for c1i = 1:length(corrAxes)
                if mxcr == 0
                    mxcr = 0.000001;
                end
                ylim(corrAxes(c1i), [-mxcr mxcr]);
                %ylim(corrAxes(c1i), [-.01 .01]);
            end
        end
        
        %PLOT LAST RM
        az = subplot(r,c,r*c-1,'Parent', parent); imagesc(az, imgaussfilt(c2.rm,1)); %clast
        title(az, sprintf('c%d',c2.ind));
        axis(az,'equal'); axis(az,'tight'); set(az,'ydir','normal','xticklabel',[],'yticklabel',[]);
        xlabel(az,sprintf('Max %.0fHz',c2.max_r));
        colormap(az,'jet');
        
        %PLOT LAST AC
        az = subplot(r,c,r*c,'Parent', parent); imagesc(az, imgaussfilt(c2.ac,1));
        hold(az,'on'); plot(az, c2.module.x, c2.module.y,'w');
        title(az, sprintf('gsc %.1f',round(c2.gridscore, 1)));
        axis(az,'equal'); axis(az,'tight'); set(az,'ydir','normal','xticklabel',[],'yticklabel',[]);
        colormap(az,'jet');
        
       
        %VERS
        gid = params.gid;
        vers = sprintf('%s_G%d_%s_deg%d_lag%.1f_blen%.2f_sig%d__win%d_',...
            v,gid,sesh,round(360/nbins),lag,binl,params.sigma,params.movmean);
        titl = sprintf('%s_moving_direction_xcorr_ncells%d',vers,length(good));
%         set(suptitle(titl),'Interpreter', 'none'); %PUT SUP TITLE AFTER ALL SUBPLOT COMMANDS << test
%         handles.fig = h;
%         handles.ax = gca;
    else
        vers = '';
        %close(h); bad??????????
    end
end

