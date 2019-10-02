function [m, h] =  driftPlot(m)
    gg = [m.g(1)];
    for i = 1:length(m.g)
        if m.g(i).before.gridscore > m.grid_thresh
            gg(end+1) = m.g(i);
        end
    end
    m.g = gg(2:end);
    t = '';
    % GENERATE FRAMES FOR FIRST TIME
    if isempty(m.windows) %~isfield(m, 'windows')
        m.windows = {}; %delete
        for i = 1:length(m.g)
            if strcmp(m.sesh,'before')
                s = m.g(i).before;
            elseif strcmp(m.sesh,'midall')
                s = m.g(i).midall;
            elseif strcmp(m.sesh,'after')
                s = m.g(i).after;
            else
                s = m.g(i).middle{m.sesh};
                t = ['mid' num2str(m.sesh)];
            end
            %   driftWindow()
            m.windows{i} = driftWindow(s, m.winsecs);
        end
        
        for u = 1:length(m.windows{1})
            %figure; imagesc(m.windows{2}{u}.rm); set(gca,'ydir','normal') 
        end
        
        
        m.crossed = cell(length(m.windows{1})+1 , length(m.g)-1, length(m.g)); % +1 is for mean frame
        for e = 1: length(m.windows{1})
            for i = 1:length(m.g)-1
                for j = i+1:length(m.g)
                    m.crossed{e,i,j} = xcorr2(m.windows{j}{e}.rm, m.windows{i}{e}.rm); %REVERSE ORDER FOR POINT OF REFERENCE 
                    %m.crossed{e,i,j} = Cross_Correlation(m.windows{i}{e}.rm, m.windows{j}{e}.rm);
                    if sum(size(m.crossed{end,i,j}) == 0) %initialize mean frame for each plot
                        m.crossed{end,i,j} = m.crossed{e,i,j};
                    else
                        %[e i j; size(m.crossed{end,i,j}) 0 ; size(m.crossed{e,i,j}) 0]
                        %mins = min( size(m.crossed{end,i,j}), size(m.crossed{e,i,j}) ); %amazingly this works
                        maxs = max( size(m.crossed{end,i,j}), size(m.crossed{e,i,j}) ); %amazingly this works
                        %running average = (previous_mean*(count-1)) + new_value)/count ;
                        m.crossed{end,i,j} = ...
                            (imresize(m.crossed{end,i,j},maxs)*(e-1) +  imresize(m.crossed{e,i,j},maxs)) /e;
      %(m.crossed{end,i,j}(1:mins(1),1:mins(2))*(e-1) +  m.crossed{e,i,j}(1:mins(1),1:mins(2)))/e;
                        
                    end
                    
                end
            end
        end
    end
    %code to show mean first
    m.wini = m.wini - 1;
    if m.wini == 0
        m.wini = length(m.windows{1}) + 1; 
    end
    %title
    if m.wini <= length(m.windows{1})
        m.title = ['Session-' num2str(m.sesh) sprintf(' [%ds - %ds]', ...
            round(m.windows{1}{m.wini}.start),...
            round(m.windows{1}{m.wini}.end))];
    else
        m.title = ['Session-' num2str(m.sesh) ' Mean Xcorr | Complete RM)'];
    end
    %plot it up
    h = plotFrame(m);
    
end

function h = plotFrame(m)
    
    if ~isfield(m,'parent')
        m.parent = figure('Position', [0, 0, 2000, 1000]);%(2*m+5)*120 , 140*r + height]); %
        set(gca,'LooseInset', get(gca,'TightInset'));
    end
    h = m.parent;
    c = length(m.g); r = c;  if r == 1; r=3; end; z = 0;
    
    %grid scatter
    a = [length(m.g)]; b = [length(m.g)];
    for i = 1:length(m.g)
        a(i) = m.g(i).before.gridscore;
        b(i) = m.g(i).midall.gridscore;
    end
    
    if ~isempty(m.g) && length(m.g) >= 2
        for i = 1:length(m.g)-1
            for j = i+1:length(m.g)
                z = z + 1;
                if strcmp(m.sesh,'before')
                    c1 = m.g(i).before; c2 = m.g(j).before; clast = m.g(r).before;
                elseif strcmp(m.sesh,'midall')
                    c1 = m.g(i).midall; c2 = m.g(j).midall; clast = m.g(r).midall;
                elseif strcmp(m.sesh,'after')
                    c1 = m.g(i).after; c2 = m.g(j).after; clast = m.g(r).after;
                else
                    %ADD CASE FOR BAD MIDDLE
                    c1 = m.g(i).middle{m.sesh}; c2 = m.g(j).middle{m.sesh};
                    clast = m.g(r).middle{m.sesh};
                end
                c1.ind = m.g(i).ind; c2.ind = m.g(j).ind; clast.ind =  m.g(r).ind;
                w1.ind = c1.ind; w2.ind = c2.ind; wlast.ind = clast.ind;
                if m.wini <= length(m.windows{1})
                    w1.rm    = m.windows{i}{m.wini}.rm; %m.wini
                    w2.rm    = m.windows{j}{m.wini}.rm;
                    wlast.rm = m.windows{r}{m.wini}.rm;
                else %plot aggragate
                    w1.rm    = c1.rm;
                    w2.rm    = c2.rm;
                    wlast.rm = clast.rm;
                end
                ax = subplot(r, c, r*r-(c*(i-1) + j-1),'Parent', m.parent); %   c*(i-1) + j-1
                cc = m.crossed{m.wini,i,j};%xcorr2(w1.rm,w2.rm);
                imagesc(ax, imgaussfilt(cc,2)); title(ax, sprintf('c%d x c%d',w1.ind,w2.ind)); %upside down?
                axis(ax,'equal'); axis(ax,'tight'); set(ax,'ydir','normal','xticklabel',[],'yticklabel',[]);
                xlabel(ax, 'xcorr');
                colormap(ax,'jet');
                hold(ax,'on'); plot(ax, size(cc,2)/2,size(cc,1)/2,'md','MarkerFaceColor','w','MarkerSize',7)
            end
            %plot rate mat
            ax = subplot(r,c,i*r,'Parent', m.parent); imagesc(ax, imgaussfilt(w1.rm,2));
            xlabel(ax,'RM');title(ax,sprintf('c%d',w1.ind));
            axis(ax,'equal'); axis(ax,'tight'); set(ax,'ydir','normal','xticklabel',[],'yticklabel',[]);   
            %assert(sum(isnan(w1.rm(:)))==0,'sum(isnan(w1.rm))==0'); colormap(ax, 'jet');
            %values now have nans
            %myColorMap = jet(256); myColorMap(1,:) = 1; colormap(ax, myColorMap);
        end
        %plot last rate mat ?
        ax = subplot(r,c,r*r,'Parent', m.parent); imagesc(ax, imgaussfilt(wlast.rm,2));
        xlabel(ax,'RM');title(ax, sprintf('c%d',wlast.ind));
        axis(ax,'equal'); axis(ax,'tight'); set(ax,'ydir','normal','xticklabel',[],'yticklabel',[]);
         %values now have nans
        %assert(sum(isnan(w1.rm(:)))==0,'sum(isnan(w1.rm))==0'); colormap(ax, 'jet');
        %myColorMap = jet(256); myColorMap(1,:) = 1; colormap(ax, myColorMap);
        % SCATTER of gridscore
        %ax = subplot(r,c,r*(r-1)+1,'Parent', m.parent,'Position',[0.05 0.09 0.35 0.45]);
        ax = subplot(r,c,1,'Parent', m.parent); 
        plot(ax, 1:length(m.g),[a; b],'o'); %legend(ax,{'before','muscimol'});
        xlabel(ax,' cell');ylabel(ax,'gridscore');title(ax, sprintf('%s','group gridscore after(red) muscimol'));
        axis(ax,'equal'); ylim(ax,[-2 2]);  %set(ax,'ydir','normal','xticklabel',[],'yticklabel',[]);
    else
        %empty group
    end
    
end

