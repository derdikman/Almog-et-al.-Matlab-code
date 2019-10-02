function CellTab(fig, tab, m)
gui.Window = fig;
gui.m = m;
gui.m.sesh = 'before';
gui.m.ind = 1;
gui.m.cgi = [];
g = gui.m.groups{1};
gui.m.c = g(1);
for i = 1:length(gui.m.groups)
    g = gui.m.groups{i};
    for j = 1:length(g)
        c = g(j);
        gui.m.cgi(c.ind,:) = [i j] ;
    end
end
s = setSesh(gui.m.ind);

%TOP UI
CellTabBox = uix.HBoxFlex( 'Parent', tab, 'Spacing', 3);
controlPanel = uix.Panel('Parent', CellTabBox,'Title', '' );
gui.viewPanel = uix.Panel('Parent', CellTabBox,'Title', 'Ready','fontweight','bold',...
    'FontSize', 11, 'TitlePosition','centerTop','ForegroundColor','b');
set( CellTabBox, 'Widths', [200 -1] );

%% VIEW PANEL
gui.viewContainer = uicontainer('Parent', gui.viewPanel,'backgroundcolor','k');
%% CONTROL PANEL
controlBoxV = uix.VBox( 'Parent', controlPanel, 'Spacing', 3);
t = [];
%% CELL
cellsBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 0 ); t = [t -1];
uicontrol('Style','text','Parent', cellsBoxV,'HorizontalAlignment', 'left','String', 'Cells:');
gui.cellsListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','Parent', cellsBoxV, ...
    'String', num2str((1:length(gui.m.cgi))'),'Value', 1,'Callback', @onCellsListSelection);
set( cellsBoxV, 'Heights', [20 -1] );
%% MID
midBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 3 ); t = [t 150];
uicontrol('Parent', midBoxV,'Style','text','HorizontalAlignment', ...
    'left', 'String','Muscimol Bin:');
gui.midListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','String', {'before';'midall'}, ...
    'Parent', midBoxV, 'Value', 1,'Callback', @onMidListSelection, 'Enable', 'on');
set( midBoxV, 'Heights', [20 -1] );

%% DELAY
delayBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 3 , 'Visible', 'on'); t = [t 60];
gui.delayLabel = uicontrol('Parent', delayBoxV,'Style','text',...
    'HorizontalAlignment', 'left','String', sprintf('Delay (s):%s',''));
delayBoxH = uix.HBox( 'Parent', delayBoxV,'Padding', 3, 'Spacing', 3 );
mnx = [-3 3]; step = 0.02; step = step/(mnx(2)-mnx(1)); %max-min
gui.delaySlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
    'Parent', delayBoxH,'Value', 1,'SliderStep',[step 5*step],'Callback', @onDelaySlide);
gui.delayEdit = uicontrol( 'Parent', delayBoxH, 'Style', 'edit', 'String', 1,...
    'Callback', @onDelayEdit);
addlistener(gui.delaySlider,'ContinuousValueChange',@(hObject, event) onDelaySliding(hObject, event));
set( delayBoxH, 'Widths', [-4 -1] );

%% FINAL HEIGHTS
set(controlBoxV, 'Heights', t );

%% CALLBACKS
% GROUP
    function onCellsListSelection( src, ~ )
        set(gui.viewPanel,'ForegroundColor',[0,0,1], 'Title', '');
        t = cellstr(get( src, 'String' ));
        gui.m.ind = str2double(t(get( src, 'Value' )));
        gui.m.sesh = 'before';
        s = setSesh(gui.m.ind);
        set(gui.viewPanel,'Title',s);
        set(gui.midListBox,'String',{'before';'midall';'after'},'Value', 1);
        if ~gui.m.g(1).after.exists
            set(gui.midListBox,'String', {'before';'midall'});
        end
        run()
    end
% MID
    function onMidListSelection( src, ~ )
        t = get( src, 'String');
        gui.m.sesh =  t{get(src,'Value')};
        s = setSesh(gui.m.ind);
        set(gui.viewPanel,'Title',s);
        run()
    end
% Delay
    function onDelaySlide( src, ~ )
        c = getSesh(1);
        t = round2r(get(src, 'Value'), 0.02);
        gui.m.delay = t/0.02; %to timesteps
        set(gui.delayEdit, 'String', t);
    end
    function onDelayEdit(src,~)
        c = getSesh(1);
        t = round2r(str2double(get(src,'String')), 0.02);
        gui.m.delay = t/0.02; %to timesteps;
        set(gui.delaySlider, 'Value', t);
    end
    function onDelaySliding(hObject,~)
        t = round2r(get(hObject,'Value'),0.02);
        set(gui.delayEdit, 'String', t);
    end

%% RUN %%
    function run()
        delete(findobj(tab,'type','axes'));
        gui.m.parent = gui.viewContainer;
        %ax = axes('Parent',gui.m.parent);
        c = gui.m.c; axi = []; s  = 7; ss = 2;
        if strcmp(gui.m.sesh,'midall')
           c.si = discretize(c.st, [-Inf, mean([c.pt(2:end); c.pt(1:end-1)]), +Inf]);
        else
           c.si = discretize(c.st, [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);
        end
        
        rows = 3;
        cols = 3;
        
        %XY        
        ax = subplot(rows,cols,1,'Parent',gui.m.parent); axi(end+1) = ax;
        plot(ax,c.px,c.py,'linewidth',1);hold(ax,'on');
        plot(ax,c.px(c.si),c.py(c.si),'w.','markersize',7);%,'linewidth',ss);
        set(ax,'Color','k');
        %center
        ax = subplot(rows,cols,9,'Parent',gui.m.parent); axi(end+1) = ax;
        plot(ax,c.px,c.py,'linewidth',1);hold(ax,'on');
        plot(ax,c.px(c.si),c.py(c.si),'w.','markersize',8);%,'linewidth',ss);
        set(ax,'Color','k');
        
        
        % RM PREV
        ax = subplot(rows,cols,3,'Parent',gui.m.parent); axi(end+1) = ax;
        imagesc(ax,c.rm);hold(ax,'on');
        
        % AC PREV
        ac = c.ac;
        ax = subplot(rows,cols,6,'Parent',gui.m.parent); axi(end+1) = ax;
        imagesc(ax,ac);hold(ax,'on');  %%%%
        ac(isnan(ac))=0;
        [k,l] = find(imregionalmax(ac));
        plot(ax,l,k,'m+','markersize',s,'linewidth',ss);
        [zmax,imax,zmin,imin]= Extrema2(ac);
        [i,j]=ind2sub(size(ac),imax); %NOAM
        %plot(ax,j,i,'mx','markersize',s,'linewidth',ss);
        title(ax,sprintf('gs %.2f',gridscore(ac,-1)),'color','m','fontsize',18);
        
        
        %RM RAW
        nbins = 50;
        mx = max(c.px);
        my = max(c.py);
        pyi = discretize(c.px, 0:mx/nbins:mx);
        pxi = discretize(c.py, 0:my/nbins:my);
        rmt = accumarray([pxi pyi], 1, [nbins nbins]);
        syi = discretize(c.px(c.si), 0:mx/nbins:mx);
        sxi = discretize(c.py(c.si), 0:my/nbins:my);
        rms = accumarray([sxi syi], 1, [nbins nbins]);
        t = rms-rmt; t(t<0)=0; rmt = rmt + t; %add in extra timestep ONLY when spikes occured faster than timestep
        rm = (rms ./ rmt*0.02);
        rm(isnan(rm)) = 0; %DO THIS?
        max(rm(:))
        %rm = rm/max(rm(:)); %normalize %REMOVE??       
        ax = subplot(rows,cols,2,'Parent',gui.m.parent); axi(end+1) = ax;
        %f = gaussian2d(nbins,2);
        %imagesc(ax, conv2(rm,f,'same'));
        imagesc(ax, rm);
  
        %AC Cross_Correlation
        ac = Cross_Correlation(rm, rm);
        ax = subplot(rows,cols,5,'Parent',gui.m.parent); axi(end+1) = ax;
        ac(isnan(ac)) = 0;
        %ac(nbins,nbins) = max(ac(:))+0.01; %for gridscore function
        aco = ac;
        ac = conv2(ac,gaussian2d(length(ac),2),'same');
        imagesc(ax,aco);hold(ax,'on');  %%%%
        [k,l] = find(imregionalmax(ac));
        %put all extrema points in dist
         ms = min(size(ac))/2; mss = ms/4; mss = round((-mss:1:mss) + ms);t = ac(mss,mss);
        [cenx ceny]=ind2sub(size(ac),find(ac==max(t(:))));
        dist = pdist2([l k],[cenx, ceny] ); [dist,ind]=sort(dist);
        if length(l) >= 7
            l=l(ind(2:7)); k=k(ind(2:7)); 
        end
        plot(ax, l,k,'m+','markersize',s,'linewidth',ss);
        viscircles(ax,[cenx ceny], dist(2)/2,'color','m');
        viscircles(ax,[cenx ceny], dist(7),'color','m');
        title(ax,sprintf('gs %.2f',gridscore2(aco,-1)),'color','m','fontsize',18);
        %AC CC SMOOTHED 
        ax = subplot(rows,cols,8,'Parent',gui.m.parent); axi(end+1) = ax;
        %ac = imgaussfilt(aco, 2,'FilterDomain','spatial'); %imshow same
        imagesc(ax,ac);hold(ax,'on');


        
        % AC2
        ac = xcorr2( rm);
        ac = ac/ max(ac(:)); %normalized ac
        ax = subplot(rows,cols,4,'Parent',gui.m.parent); axi(end+1) = ax;
        acg = imgaussfilt(ac, 2,'FilterDomain','spatial'); %imshow
        imagesc(ax,ac); hold(ax,'on'); %%%%
        [k,l] = find(imregionalmax(acg));
        dist = pdist2([l k],[nbins, nbins] ); [dist,ind]=sort(dist);
        if length(l) >= 7
            l=l(ind(2:7)); k=k(ind(2:7)); 
        end
        plot(ax,l,k,'m+','markersize',s,'linewidth',ss);
        %plot(ax,[nbins l(2:3)' nbins],[nbins k(2:3)' nbins],'m'); triangle
        if length(l) >= 2
            [~,ceny]=max(max(ac));[~,cenx]=max(max(ac'));
            viscircles(ax,[cenx ceny], dist(2)/2,'color','m');
            viscircles(ax,[cenx ceny], dist(2)*3/2,'color','m');
        end
        title(ax,sprintf('gs %.2f',gridscore2(ac,-1)),'color','m','fontsize',18);
        % AC2 SMOOTHED
        ax = subplot(rows,cols,7,'Parent',gui.m.parent) ; axi(end+1) = ax;
        %[cs, rs] = imfindcircles(ac>std(ac(:)),[10 40]);
        imagesc(ax,acg); %hold on;
        
               
        %imagesc(ax,imgaussfilt(rm, 1)); %imshow
        % CORRELATION
        %imagesc(corrcoef(rm));
        %imagesc(xcorr2(rm,rm));
        %a(a > 0.1) = 0.1; 
        %imagesc(a);
        
        for i = 1:length(axi)
            colormap(axi(i), 'jet');
            set(axi(i),'YDir','normal');
            axis(axi(i),'tight');
        end
        
    end %end run()

%% UTIL
    function s = setSesh(ind)
        gui.m.g = gui.m.groups{gui.m.cgi(ind,1)};
        i = gui.m.cgi(gui.m.ind,2);
        if strcmp(gui.m.sesh,'before')
            gui.m.c = gui.m.g(i).before; s = 'before';
        elseif strcmp(gui.m.sesh,'midall')
            gui.m.c = gui.m.g(i).midall; s = 'midall';
        elseif strcmp(gui.m.sesh,'after')
            gui.m.c = gui.m.g(i).after; s = 'after';
        else
            gui.m.c = gui.m.g(i).middle{gui.m.sesh}; s = ['mid' num2str(gui.m.sesh)];
        end
        s = sprintf('Cell %d %s',gui.m.ind, s);
    end
    
    function train
        train = [];%zeros(length(m.g), length(c.pt));
        c = [];
        for i = 1:length(gui.m.g)
            c = getSesh(i);
            if strcmp(gui.m.sesh,'before')
                edges = [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf];
            elseif strcmp(gui.m.sesh,'midall')
                edges = [-Inf, mean([c.pt(2:end); c.pt(1:end-1)]), +Inf];
            elseif strcmp(gui.m.sesh,'after')
                edges = [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf];
            else
                c = gui.m.g(i).middle{gui.m.sesh};
                edges = [-Inf, mean([c.pt(2:end); c.pt(1:end-1)]), +Inf];
            end
            %if isnumeric(gui.m.sesh)
            if isnan(edges(2))
                edges(2) = 1.1;
            end
            I = discretize(c.st, edges);
            train(i,:) = zeros(1, length(c.pt));
            train(i,I) = 1;
        end
        gui.m.train = logical(train);
    end

run();
end

function f = gaussian2d(N,sigma)
  % N is grid size, sigma speaks for itself
 [x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
 f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
 f=f./sum(f(:));
end














