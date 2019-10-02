function TabCell(fig, tab, m)
v.Window = fig;
m.sesh = 'before';
m.ind = 1;
m.cgi = [];
g = m.groups{1};
m.c = g(1);
for i = 1:length(m.groups)
    g = m.groups{i};
    for j = 1:length(g)
        c = g(j);
        m.cgi(c.ind,:) = [i j] ;
    end
end
s = setSesh();
m.t0 = 1;
m.tn = length(m.c.pt);
m.f = 500;
m.win = 0;

%#ok<*NASGU>

%TOP UI
CellTabBox = uix.HBoxFlex( 'Parent', tab, 'Spacing', 3);
controlPanel = uix.Panel('Parent', CellTabBox,'Title', '' );
v.viewPanel = uix.Panel('Parent', CellTabBox,'Title', 'Ready','fontweight','bold',...
    'FontSize', 11, 'TitlePosition','centerTop','ForegroundColor','b');
set( CellTabBox, 'Widths', [200 -1] );

%% VIEW PANEL
v.viewContainer = uicontainer('Parent', v.viewPanel,'backgroundcolor','k');
%% CONTROL PANEL
controlBoxV = uix.VBox( 'Parent', controlPanel, 'Spacing', 3);
t = []; pad = 1; spc = 1; sld = 35;
%% t0
t0BoxV = uix.VBox( 'Parent', controlBoxV,'Padding',pad,'Spacing',spc, 'Visible', 'on'); t = [t sld];
v.t0Label = uicontrol('Parent', t0BoxV,'Style','text',...
    'HorizontalAlignment', 'left','String', sprintf('t0 (s):%s',''));
t0BoxH = uix.HBox( 'Parent', t0BoxV,'Padding',pad,'Spacing',spc );
mnx = [1 floor(length(m.c.pt)*0.02)]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
v.t0Slider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
    'Parent', t0BoxH,'Value', 1,'SliderStep',[step 10*step],'Callback', @ont0Slide);
v.t0Edit = uicontrol( 'Parent', t0BoxH, 'Style', 'edit', 'String', 1,'Callback', @ont0Edit);
addlistener(v.t0Slider,'ContinuousValueChange',@(hObject, event) ont0Sliding(hObject, event));
set( t0BoxH, 'Widths', [-4 -1] );
%% tn
tnBoxV = uix.VBox( 'Parent', controlBoxV,'Padding',pad,'Spacing',spc, 'Visible', 'on'); t = [t sld];
v.tnLabel = uicontrol('Parent', tnBoxV,'Style','text',...
    'HorizontalAlignment', 'left','String', sprintf('tn (s):%s',''));
tnBoxH = uix.HBox( 'Parent', tnBoxV,'Padding',pad,'Spacing',spc );
mnx = [1 floor(length(m.c.pt)*0.02)]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
v.tnSlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
    'Parent', tnBoxH,'Value', mnx(2), 'SliderStep',[step 10*step],'Callback', @ontnSlide);
v.tnEdit = uicontrol( 'Parent', tnBoxH, 'Style', 'edit', 'String', mnx(2),'Callback', @ontnEdit);
addlistener(v.tnSlider,'ContinuousValueChange',@(hObject, event) ontnSliding(hObject, event));
set( tnBoxH, 'Widths', [-4 -1] );
%% f
fBoxV = uix.VBox( 'Parent', controlBoxV,'Padding',pad,'Spacing',spc, 'Visible', 'on'); t = [t sld];
v.fLabel = uicontrol('Parent', fBoxV,'Style','text',...
    'HorizontalAlignment', 'left','String', sprintf('f (s):%s',''));
fBoxH = uix.HBox( 'Parent', fBoxV,'Padding',pad,'Spacing',spc );
mnx = [1 500]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
v.fSlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
    'Parent', fBoxH,'Value', m.f,'SliderStep',[step 10*step],'Callback', @onfSlide);
v.fEdit = uicontrol( 'Parent', fBoxH, 'Style', 'edit', 'String', m.f,'Callback', @onfEdit);
addlistener(v.fSlider,'ContinuousValueChange',@(hObject, event) onfSliding(hObject, event));
set( fBoxH, 'Widths', [-4 -1] );
%% win
winBoxV = uix.VBox( 'Parent', controlBoxV,'Padding',pad,'Spacing',spc, 'Visible', 'on'); t = [t sld];
v.winLabel = uicontrol('Parent', winBoxV,'Style','text',...
    'HorizontalAlignment', 'left','String', sprintf('win (s):%s',''));
winBoxH = uix.HBox( 'Parent', winBoxV,'Padding',pad,'Spacing',spc );
mnx = [0 1e4]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
v.winSlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
    'Parent', winBoxH,'Value', m.f, 'SliderStep',[step 10*step],'Callback', @onwinSlide);
v.winEdit = uicontrol( 'Parent', winBoxH, 'Style', 'edit', 'String', m.f,'Callback', @onwinEdit);
addlistener(v.winSlider,'ContinuousValueChange',@(hObject, event) onwinSliding(hObject, event));
set( winBoxH, 'Widths', [-4 -1] );
%% BUTTONS
    buttonBoxH = uix.HBox( 'Parent', controlBoxV, 'Spacing', spc); t = [t 30];
    uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'gs>thr','Callback', @gsb );
    uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'hd>0.4','Callback', @hdb,'Enable','On');
    uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'non','Callback', @nonb,'Enable','On');
    uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'gs><','Callback', @decb,'Enable','On');
    uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'all','Callback', @allb,'Enable','On');
%% MID
midBoxV = uix.VBox( 'Parent', controlBoxV,'Padding',pad,'Spacing',spc); t = [t 70];
uicontrol('Parent', midBoxV,'Style','text','HorizontalAlignment', ...
    'left', 'String','Muscimol Bin:');
midListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','String', {'before';'midall'}, ...
    'Parent', midBoxV, 'Value', 1,'Callback', @onMidListSelection, 'Enable', 'on');
set( midBoxV, 'Heights', [20 -1] );
%% CELL
cellsBoxV = uix.VBox( 'Parent', controlBoxV,'Padding',pad,'Spacing',spc); t = [t -1];
uicontrol('Style','text','Parent', cellsBoxV,'HorizontalAlignment', 'left','String', 'Cells:');
v.cellsListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','Parent', cellsBoxV, ...
    'String', num2str((1:length(m.cgi))'),'Value', 1,'Callback', @onCellsListSelection);
set( cellsBoxV, 'Heights', [20 -1] );

%% FINAL HEIGHTS
set(controlBoxV, 'Heights', t );

%% CALLBACKS
function gsb(~, ~)
    t = [];
    for i = 1:length(m.cgi)
        g = m.groups{m.cgi(i,1)};
        k = m.cgi(i,2);
        if  strcmp(m.sesh,'before')
            c = g(k).before;
        else
            c = g(k).midall;
        end
        if c.gridscore > v.Window.UserData.gridThresh
            t = [t i];
        end
    end
    set(v.cellsListBox, 'String',t); 
end
function hdb(~, ~)
        t = [];
    for i = 1:length(m.cgi)
        g = m.groups{m.cgi(i,1)};
        k = m.cgi(i,2);
        if  strcmp(m.sesh,'before')
            c = g(k).before;
        else
            c = g(k).midall;
        end
        if c.rayleigh_score > 0.7
            t = [t i];
        end
    end
    set(v.cellsListBox, 'String',t);
end
function nonb(~, ~)
        t = [];
        v.Window.UserData.gridThreshMid
        m.sesh
    for i = 1:length(m.cgi)
        g = m.groups{m.cgi(i,1)};
        k = m.cgi(i,2);
        if  strcmp(m.sesh,'before')
            c = g(k).before;
        else
            c = g(k).midall;
        end
        if c.gridscore < v.Window.UserData.gridThreshMid...
                && c.rayleigh_score < 0.3
            t = [t i];
        end
    end
    set(v.cellsListBox, 'String',t);
end
function decb(~, ~)
     t = [];
    for i = 1:length(m.cgi)
        g = m.groups{m.cgi(i,1)};
        k = m.cgi(i,2);
        if g(k).midall.gridscore < v.Window.UserData.gridThreshMid...
        && g(k).before.gridscore > v.Window.UserData.gridThresh ...
        && len(g(k).midall.st) > 700  && len(g(k).before.st) > 700
            t = [t i];
        end
    end
    
    set(v.cellsListBox, 'String',t);
end
function allb(~, ~)
    set(v.cellsListBox, 'String',[1:301]);
end
% CELL TAB
    function onCellsListSelection( src, ~ )
        set(v.viewPanel,'ForegroundColor',[0,0,1], 'Title', '');
        t = cellstr(get( src, 'String' ));
        m.ind = str2double(t(get( src, 'Value' )));
        %m.sesh = 'before'; %PUT BACK IN
        s = setSesh();
        set(v.viewPanel,'Title',s);
        set(midListBox,'String',{'before';'midall';'after'},'Value', 1);
        if ~m.g(1).after.exists
            set(midListBox,'String', {'before';'midall'});
        end
        if m.tn > length(m.c.pt)
            m.tn = length(m.c.pt);
        end
        m.t0 = 1;
        m.tn = length(m.c.pt);
        mnx = [1 floor(length(m.c.pt)*0.02)]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
        set(v.tnSlider,'Min', mnx(1),'Max', mnx(2),'Parent', tnBoxH,'Value',  floor(length(m.c.pt)*0.02) );
        run()
    end
% MID
    function onMidListSelection( src, ~ )
        t = get( src, 'String');
        m.sesh = t{get(src,'Value')};
        s = setSesh();
        m.t0 = 1;
        m.tn = length(m.c.pt);
        mnx = [1 floor(length(m.c.pt)*0.02)]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
        set(v.tnSlider,'Min', mnx(1),'Max', mnx(2),'Parent', tnBoxH,'Value', floor(length(m.c.pt)*0.02) );
        set(v.viewPanel,'Title',s);
        run()
    end
% t0
    function ont0Slide( src, ~ )
        t = round2r(get(src, 'Value'), 1);
        m.t0 = t/0.02; %to timesteps
        set(v.t0Edit, 'String', t);
        run();
    end
    function ont0Edit(src,~)
        t = round2r(str2double(get(src,'String')), 1);
        m.t0 = t/0.02; %to timesteps;
        set(v.t0Slider, 'Value', t);
        run();
    end
    function ont0Sliding(hObject,~)
        t = round2r(get(hObject,'Value'),1);
        set(v.t0Edit, 'String', t);
    end
% tn
    function ontnSlide( src, ~ )
        t = round2r(get(src, 'Value'), 1);
        m.tn = t/0.02; %to timesteps
        set(v.tnEdit, 'String', t);
        run();
    end
    function ontnEdit(src,~)
        t = round2r(str2double(get(src,'String')), 1);
        m.tn = t/0.02; %to timesteps;
        set(v.tnSlider, 'Value', t);
        run();
    end
    function ontnSliding(hObject,~)
        t = round2r(get(hObject,'Value'),1);
        set(v.tnEdit, 'String', t);
    end
% f
    function onfSlide( src, ~ )
        t = round2r(get(src, 'Value'), 1);
        m.f = t; %to timesteps
        set(v.fEdit, 'String', t);
        run();
    end
    function onfEdit(src,~)
        t = round2r(str2double(get(src,'String')), 1);
        m.f = t; %to timesteps;
        set(v.fSlider, 'Value', t);
        run();
    end
    function onfSliding(hObject,~)
        t = round2r(get(hObject,'Value'),1);
        set(v.fEdit, 'String', t);
    end
% win
    function onwinSlide( src, ~ )
        t = round2r(get(src, 'Value'), 1);
        m.win = t; %to timesteps
        set(v.winEdit, 'String', t);
        run();
    end
    function onwinEdit(src,~)
        t = round2r(str2double(get(src,'String')), 1);
        m.win = t; %to timesteps;
        set(v.winSlider, 'Value', t);
        run();
    end
    function onwinSliding(hObject,~)
        t = round2r(get(hObject,'Value'),1);
        set(v.winEdit, 'String', t);
    end



%% UTIL
    function s = setSesh()
        m.g = m.groups{m.cgi(m.ind,1)};
        i = m.cgi(m.ind,2);
        if strcmp(m.sesh,'before')
            m.c = m.g(i).before; s = 'before';
        elseif strcmp(m.sesh,'midall')
            m.c = m.g(i).midall; s = 'midall';
        elseif strcmp(m.sesh,'after')
            m.c = m.g(i).after; s = 'after';
        else
            m.c = m.g(i).middle{m.sesh}; s = ['mid' num2str(m.sesh)];
        end
        s = sprintf('Cell %d g%di%d %s',m.ind,m.cgi(m.ind,1),m.cgi(m.ind,2),s);
    end


    function f = gaussian2d(N,sigma)
        % N is grid size, sigma speaks for itself
        [x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
        f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
        f=f./sum(f(:));
    end

%% RUN %%
    function run()
        delete(findobj(tab,'type','axes'));
        tfs = 12; ms  = 7; lw = 2;
        rows = 3; cols = 3; axi = [];
        m.parent = v.viewContainer;
        if m.t0 > m.tn
            m.tn = tui.m.t0 + 1;
        end
        %c.ind
        nbins = 200;
        
        c = m.c; dt = round(median(diff(c.pt)),3); %to ms
        %a = patchTrajectoryLinear(c.pt,c.px,c.py,dt,dt*1.9); %NO NEED NOW?
        pw = m.t0:m.tn; %pw = 1:length(c.pt); %DO PW HERE
        c.si = discretize(c.st, [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);
        sw = c.si(c.si>m.t0 & c.si< m.tn);
        mx = max(c.px(pw));my = max(c.py(pw));
        pxi = discretize(c.px(pw), 0:mx/nbins:mx);
        pyi = discretize(c.py(pw), 0:my/nbins:my);
        %t = diff(c.pt(pw)); t = [median(t); t];
        rmt = accumarray([pyi pxi], 1, [nbins nbins]);
        rmt(rmt<1) = 1e10;
        %rmt(rmt<min(t)) = 1e10;
        sxi = discretize(c.px(sw), 0:mx/nbins:mx);
        syi = discretize(c.py(sw), 0:my/nbins:my);
        rms = accumarray([syi sxi], 1, [nbins nbins]);
        rm = rms./(rmt*dt); [mm, i] = max(rm(:));
        rm(isnan(rm)) = 0; %DO THIS? %ind2sub(size(rm),i);
        %t = rms-rmt; t(t<0)=0; rmt = rmt + t; %add in extra timestep ONLY when spikes occured faster than timestep
        %SMOOTHED TRAIN
        tr = createMsSpikeTrain(c.st(c.st>(c.pt(m.t0)) & c.st< c.pt(m.tn)));
        trs = tr;
        if m.win > 0
            trs = movmean(trs,m.win);
        end
        st = [1:len(trs)]/1000;
        sti = discretize(st, [-Inf,mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);%puts time into c.pt units
        sxsi = discretize(c.px(sti), 0:mx/nbins:mx);
        sysi = discretize(c.py(sti), 0:my/nbins:my);
        rmss = accumarray([sysi sxsi], trs, [nbins nbins]);
        %         pxi = discretize(c.px, 0:mx/nbins:mx);
        %         pyi = discretize(c.py, 0:my/nbins:my);
        %         %t = diff(c.pt(pw)); t = [median(t); t];
        %         rmt = accumarray([pyi pxi], 1, [nbins nbins]);
        rm2 = rmss./rmt;
        
        %%%%%%%%%%
        %PLOTS
        %%%%%%%%%%
        
        show=[1 1 1 1 1 1 1 1 1];
        
        %PLOT 1
        if show(1)
        ax = subplot(rows,cols,1,'Parent',m.parent); axi(end+1) = ax;
        %XY
        plot(ax,c.px(pw),c.py(pw),'b','linewidth',0.3);hold(ax,'on');%set(ax,'color','w');
        plot(ax,c.px(sw),c.py(sw),'w.','markersize',3);%,'linewidth',ss);
        title(ax,sprintf('spikes=%d',length(sw)),'color','m','fontsize',tfs);
        end
        
        %PLOT 2
        if show(2)
        ggg = 7;
        ax = subplot(rows,cols,2,'Parent',m.parent);axi(end+1) = ax;
        ac2 = xcorr2(rm2);
        imagesc(ax,imgaussfilt(ac2,ggg));colormap jet; axis(gca,'square');
        %plot(trs)
        set(gca,'YDir','normal'); axis(gca,'tight');
        title(sprintf('gs %.2f',gridscore2(ac2,ggg)),'color','m','fontsize', tfs);
        end
        
        %PLOT 3
        if show(3)
        ax = subplot(rows,cols,3,'Parent',m.parent);axi(end+1) = ax;
        imagesc(ax,imgaussfilt(rm2,1));colormap jet; axis(gca,'square');
        %plot(trs)
        set(gca,'YDir','normal'); axis(gca,'tight');
        title(sprintf('win %.dms',m.win),'color','m','fontsize', tfs);
        end
        
        % PLOT 4
        if show(4)
        ax = subplot(rows,cols,4,'Parent',m.parent); axi(end+1) = ax; cla(ax);
        %FOURIER
        %          spkms = floor(c.st*1000)+1; %to millisecs
        %          spkms = spkms-min(spkms)+100; %add 100 leading zeros
        %          train= zeros(1,max(spkms)+100);%add 100 trailing zeros
        %          train(spkms)=1;
        %t = zeros(c.si(end)+10,1); %add 10 trailing zeros
        %t(c.si)=1;
        plotfft(ax,tr,1000,m.f,100); %time in ms;
        hold(ax, 'on');
        plotfft(ax,trs,1000,m.f,100);
        title(ax,'FOURIER');
        end
        
        %PLOT 5
        if show(5)
        ax = subplot(rows,cols,5,'Parent',m.parent); axi(end+1) = ax;
        %RM
        %maxfiring rate
        %rm = imgaussfilt(rm,1);
        %rm = rm/max(rm(:)); %normalize %REMOVE??
        %f = gaussian2d(nbins,2);
        %imagesc(ax, conv2(rm,f,'same'));
        imagesc(ax, imgaussfilt(rm,1)); %rm
        title(ax,sprintf('max %.1fHz',mm),'color','m','fontsize', tfs);
        end
            

        
        %PLOT 6
        if show(6)
        ax = subplot(rows,cols,6,'Parent',m.parent); axi(end+1) = ax;
        % AC
        ac = xcorr2(rm);
        ac = ac/ max(ac(:)); %normalized ac;
        acg = imgaussfilt(ac, 2,'FilterDomain','spatial'); %imshow
        imagesc(ax,acg); hold(ax,'on'); %%%%
        [k,l] = find(imregionalmax(acg));
        dist = pdist2([l k],[nbins, nbins] ); [dist,ind]=sort(dist);
        if length(l) >= 7
            l=l(ind(2:7)); k=k(ind(2:7));
        end
        %plot(ax,l,k,'m+','markersize',s,'linewidth',ss);
        acg = imgaussfilt(ac, 3,'FilterDomain','spatial'); %imshow
        [k,l] = find(imregionalmax(acg));
        dist = pdist2([l k],[nbins, nbins] ); [dist,ind]=sort(dist);
        if length(l) >= 7
            l=l(ind(2:7)); k=k(ind(2:7));
        end
        plot(ax,l,k,'yx','markersize',ms,'linewidth',lw);
        %plot(ax,[nbins l(2:3)' nbins],[nbins k(2:3)' nbins],'m'); triangle
        if length(l) >= 2
            [~,ceny]=max(max(ac));[~,cenx]=max(max(ac'));
            %viscircles(ax,[cenx ceny], dist(2)/2,'color','m'); %ADD BACK
            %viscircles(ax,[cenx ceny], dist(2)*3/2,'color','m');
        end
        %rms(rms>0)=1;
        title(ax,sprintf('gs%.2f gsa%.2f gso%.2f',gridscore2(ac,2),...
            gridscore2(xcorr2(imgaussfilt(rm,1)),2),...
            c.gridscore),'color','m','fontsize',tfs); %rms
        %         acg = imgaussfilt(ac, 3,'FilterDomain','spatial');
        %         %acg = imresize(acg, [70 120]);
        %         imagesc(ax,acg); hold(ax,'on'); %%%%
        %         cent = size(acg)/2;
        %         plot(ax,cent(2),cent(1),'wx','markersize',10);
        module = Find_Module(acg);
        plot(ax,module.hex_peaks(:,1),module.hex_peaks(:,2),'ro');
        %plot(ax,module.hex_peaks2(:,1),module.hex_peaks2(:,2),'y*');
        plot(ax,module.x,module.y,'r');
        %title(ax,sprintf(' %s','Module'),'color','m','fontsize',tfs);
        % ax = subplot(rows,cols,3,'Parent',m.parent)  %PLOT 1 TITLE
        %title(ax,sprintf('gs %.2f \ngs2 %.2f \ngs3 %.2f',gridscore2(ac,2),...
        % AC2 SMOOTHED
        %[cs, rs] = imfindcircles(ac>std(ac(:)),[10 40]);
        end
        
        % PLOT 7
        if show(7)
        ax = subplot(rows,cols,7,'Parent',m.parent); axi(end+1) = ax;
        % TRAIN
        t = zeros(c.si(end),1);
        t(c.si)=1;
        %          bin_s = 0.1;
        %          t = histcounts(c.st,0:bin_s:c.st(end));
        %          plot(ax,t,'y');
        
        tac = xcorr(tr);%xcorr(decimate(tr,1e3)); st = round(len(tac)/2);
        %tac = tac(st:st+1e5);
        %tac = movmean(tac,50);
        %tac(1:1000) = 0;
        plot(ax,tac(round(len(tac)/2)+5e2:end));
        hold(ax, 'on');
        %tac = tac(st:st+1e5);
        tac = xcorr(trs);%xcorr(decimate(trs,1e3));
        plot(ax,tac(round(len(tac)/2)+5e2:end));
        %plot(ax,tac);
        %title(ax,sprintf('Train bin=%.2fs', 1),);
        title(ax,'time acorr','color','m','fontsize',tfs);
        end
        
        % PLOT 8
        if show(8)
        ax = subplot(rows,cols,8,'Parent',m.parent); axi(end+1) = ax;
        % RAYLEIGH HD
%         d = m.c;
%         rbin = linspace(-pi,pi,1000)
%         rd = histcounts(rad2deg(c.hd(c.si))+180,0:3:360)./histcounts(rad2deg(c.hd)+180,0:3:360);
%         rd = smooth(rd,15);
        
        %         sti = discretize(st, [-Inf,mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);%puts time into c.pt units
        %         sxsi = discretize(c.px(sti), 0:mx/nbins:mx);
        %         sysi = discretize(c.py(sti), 0:my/nbins:my);
        %         rmss = accumarray([sysi sxsi], trs, [nbins nbins]);
        
        rbin = linspace(-pi,pi,100);
        shdi = discretize(wrapToPi(c.hd(sti)), rbin);
        rd2 = accumarray(shdi, trs,[len(rbin) 1])';
        rd = accumarray(discretize(wrapToPi(c.hd(c.si)), rbin), 1,[len(rbin) 1])';
        rdt = accumarray(discretize(wrapToPi(c.hd), rbin), 1, [len(rbin) 1])';
        %rbin = -2*pi:0.1:2*pi;
        %shdi = discretize(c.hd(sti),rbin);
        %rd2 = accumarray(shdi, trs,[length(rbin) 1]);
        rs = m.c.rayleigh_score;
        rd2a = rd2./rdt;  rda = rd./rdt;
        rd2a(isnan(rd2a)) = 0; rda(isnan(rda)) = 0;
        
        r2 = sqrt( sum(cos(rbin).*rd2a).^2+sum(sin(rbin).*rd2a).^2)/sum(rd2a);
        %norm_val = sqrt(sin(c.hd).^2+cos(c.hd).^2) %nansum(c.hd); %denominator to normalize to 1
        %veclen = sqrt(sum(xbin/len(xbin))^2+sum(ybin/len(ybin))^2); %length of x and y component of each bin
        %rs2 = veclen/sqrt(sum(cos(c.hd)/len(c.hd))^2+sum(sin(c.hd)/len(c.hd))^2);rayleigh_score(c.hd,c.hd(c.si)))
        %plot(ax, (rad2deg(rbin(2:end))),max(rd)*rd2/max(rd2),'g','linewidth',0.7); xlim(ax,[0 360]);
        plot(ax, rbin,rda,'g','linewidth',0.7);
        hold on;
        plot(ax, rbin,rd2a,'y','linewidth',0.7);
        title(ax,sprintf('HD ray orig=%.2f smth=%.2f',c.rayleigh_score,r2),'color','m','fontsize',tfs);
        end
        
       
        
        %PLOT 9
        if show(9)
        ax = subplot(rows,cols,9,'Parent',m.parent); %axi(end+1) = ax;
        % ISI
        bins = 10.^[-4:0.05:2];
  %      histogram(ax,diff(m.c.st),bins); 
        plot(ax,histcounts(diff(m.c.st),bins),'r','linewidth',lw);%imagesc(ax, rm); %rm
        ax.XTick = round(linspace(1,length(bins),10));
        %l = cellstr(num2str(bins(ax.XTick')));
        ax.XTickLabel= round(bins(ax.XTick),4);
        ax.XTickLabelRotation = 270;
        title(ax,sprintf('ISI (s)',mm),'color','m','fontsize', tfs);
        set(ax,'Ycolor','m');
        set(ax,'Xcolor','m');
        set(ax,'Color','k');
        %axis(ax,'square');
        end
        
        
        % PLOT 9
        %ax = subplot(rows,cols,9,'Parent',m.parent); axi(end+1) = ax; cla(ax);
        %FOURIER
        %          spkms = floor(c.st*1000)+1; %to millisecs
        %          spkms = spkms-min(spkms)+100; %add 100 leading zeros
        %          train= zeros(1,max(spkms)+100);%add 100 trailing zeros
        %          train(spkms)=1;
        %t = zeros(c.si(end)+10,1); %add 10 trailing zeros
        %t(c.si)=1;
        
        %plotfft(ax,fft(ftrain),1000); %time in ms;
        %title(ax,'FOURIER BANDPASS');
        
        %          % PLOT 9
        %          ax = subplot(rows,cols,9,'Parent',m.parent); axi(end+1) = ax;
        %          % TIME DIFF
        %          plot(ax, c.pt);
        %          hold(ax,'on'); plot(ax,c.si,c.st,'r.');
        %          df = diff(c.pt);
        %          hold(ax,'off'); hist(ax,df,10);
        %          %title(ax,sprintf('Time diff mx%.2f m%.3f s%.2f', max(df), mean(df), std(df) ),'color','m','fontsize',tfs-2);
        %          title(ax,sprintf('%s','Time diff (s)'),'color','m','fontsize',tfs);
        
        % FOR PAPER - REMOVE
        %          ax = subplot(rows,cols,2,'Parent',m.parent); axi(end+1) = ax; cla(ax);
        %          % AC Clean
        %          acg = imgaussfilt(ac, 3,'FilterDomain','spatial');imagesc(ax,acg);
        %          title(ax,sprintf('gridscore=%.1f',round(gridscore2(ac,3),1)),'color','m','fontsize',tfs); %rms
        %
        %          ax = subplot(rows,cols,3,'Parent',m.parent); axi(end+1) = ax; cla(ax);
        %          % RAYLEIGH
        %          d = m.c;
        %          rd = histcounts(rad2deg(c.hd(c.si))+180,0:3:360)./histcounts(rad2deg(c.hd)+180,0:3:360);
        %          rd = smooth(rd,45);
        %          plot(ax, 3:3:360,rd,'g','linewidth',3);
        %          rs = d.rayleigh_score;
        %          title(ax,sprintf('Rayleigh=%.1f',round(rs,1)),'color','m','fontsize',tfs);
        %          xlabel(ax,'angle','color','m','fontsize',tfs-2)
        
        % PLOT 3
        %{
         ax = subplot(rows,cols,3,'Parent',m.parent);axi(end+1) = ax;
         sig = 1.4;
         ac = xcorr2( rm);
         ac = ac/ max(ac(:)); %normalized ac;
         acg = imgaussfilt(ac, sig,'FilterDomain','spatial'); %imshow
         imagesc(ax,acg); hold(ax,'on'); %%%%
         module = Find_Module(acg);
         plot(ax,module.hex_peaks(:,1),module.hex_peaks(:,2),'ro');
        %plot(ax,module.hex_peaks2(:,1),module.hex_peaks2(:,2),'y*');
         plot(ax,module.x,module.y,'r');
         title(ax,sprintf('gs %.2f', gridscore2(ac,sig)),'color','m','fontsize',tfs);
%         % reconstruct        %floor(c.st*1000)+1;
%         newst = 1:length(ftrain); %newst(~train)=0;
%         newst = newst(boolean( round( real(ftrain)) ))/1000;
%         newsi = discretize(newst, [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf]);
        %sw = c.si(c.si>m.t0 & c.si< m.tn);
        %[newrm, ~] = Create_Rate_Map(c.px, c.py, c.pt, c.px(newsi), c.py(newsi), newst, true, 50);
        %plot(ax,c.px(newsi),c.py(newsi),'w.','markersize',3);%,'linewidth',ss);
        %title(ax,sprintf('spikes=%d',length(sw)),'color','m','fontsize',tfs);
        %}
        
        
        
        %ALLLL
        for i = 1:length(axi)
            colormap(axi(i), 'jet');
            set(axi(i),'YDir','normal');
            axis(axi(i),'tight');
            set(axi(i),'Ycolor','m');
            set(axi(i),'Xcolor','m');
            set(axi(i),'Color','k');
            axis(axi(i),'square');
        end
        
    end %end run()


end



%          % PLOT 1
%          ax = subplot(rows,cols,1,'Parent',m.parent); axi(end+1) = ax; %cla(ax);
%          %FOURIER
%          spkms = floor(c.st*1000)+1; %to millisecs
%          %spkms = spkms-min(spkms)+100; %add 100 leading zeros
%          %train= zeros(1,max(spkms)+100);%add 100 trailing zeros
%          train= zeros(1,length(spkms));
%          train(spkms)=1;          %t = zeros(c.si(end)+10,1); %add 10 trailing zeros        %t(c.si)=1;
%          %plotfft(ax,train,1000); %time in ms;
% %          y=fft(train-mean(train));y(1:500)=0; %deletes first few;
% %          l = length(y); l2 = floor(l/2);
% %          y = abs(y(1:l2));
% %          r = 500;y = decimate(y,r); %y=abs(y);
%          f1 = 1e-9; f2 = 500-1e-9; fs = 1000;
%          ftrain = train;
%          %ftrain = bandpassfft(train, f1, f2, fs);
%          y=fft(train-mean(ftrain));y(1:500)=0; %deletes first few;
%          l = length(y); l2 = floor(l/2);
%          y = abs(y(1:l2));
%          f = 1000*(1:length(y))/(length(y)*2);
%         % y = smooth(y,1e2);%'sgolay',
%          fi = round((length(y)*    20    *2/1000));
%          plot(ax,f(1:fi),smooth(y(1:fi),500),'linewidth',3); ax.YLim = [0 100]; hold(ax,'on');
%          plot(ax,[7.9999 8],[ylim(ax)],'r.','linewidth',1.2,'markersize',15);  %max(y(1:round(l/2)) )
%          xlabel(ax,'Hz');
%          title(ax,sprintf('gs %.1f  len %d ',round(c.gridscore,1),length(y)),'color','m','fontsize',tfs); ax.YLim = [0 100]; %hold(ax,'off');
%          %figure(2);
%         % plotfft(gca,train,1000);






%center
%ax = subplot(rows,cols,9,'Parent',m.parent); axi(end+1) = ax;
%         plot(ax,c.px,c.py,'linewidth',1);hold(ax,'on');
%         plot(ax,c.px(c.si),c.py(c.si),'w.','markersize',8);%,'linewidth',ss);
%         set(ax,'Color','k');


%         % RM PREV
%         ax = subplot(rows,cols,3,'Parent',m.parent); axi(end+1) = ax;
%         imagesc(ax,c.rm);hold(ax,'on');
%
%         % AC PREV
%         ac = c.ac;
%         ax = subplot(rows,cols,6,'Parent',m.parent); axi(end+1) = ax;
%         imagesc(ax,ac);hold(ax,'on');  %%%%
%         ac(isnan(ac))=0;
%         [k,l] = find(imregionalmax(ac));
%         plot(ax,l,k,'m+','markersize',s,'linewidth',ss);
%         [zmax,imax,zmin,imin]= Extrema2(ac);
%         [i,j]=ind2sub(size(ac),imax); %NOAM
%         %plot(ax,j,i,'mx','markersize',s,'linewidth',ss);
%         title(ax,sprintf('gs %.2f',gridscore(ac,-1)),'color','m','fontsize',18);




%         %AC Cross_Correlation
%         ac = Cross_Correlation(rm, rm);
%         ax = subplot(rows,cols,5,'Parent',m.parent); axi(end+1) = ax;
%         ac(isnan(ac)) = 0;
%         %ac(nbins,nbins) = max(ac(:))+0.01; %for gridscore function
%         aco = ac;
%         ac = conv2(ac,gaussian2d(length(ac),2),'same');
%         imagesc(ax,aco);hold(ax,'on');  %%%%
%         [k,l] = find(imregionalmax(ac));
%         %put all extrema points in dist
%          ms = min(size(ac))/2; mss = ms/4; mss = round((-mss:1:mss) + ms);t = ac(mss,mss);
%         [cenx ceny]=ind2sub(size(ac),find(ac==max(t(:))));
%         dist = pdist2([l k],[cenx, ceny] ); [dist,ind]=sort(dist);
%         if length(l) >= 7
%             l=l(ind(2:7)); k=k(ind(2:7));
%         end
%         plot(ax, l,k,'m+','markersize',s,'linewidth',ss);
%         viscircles(ax,[cenx ceny], dist(2)/2,'color','m');
%         viscircles(ax,[cenx ceny], dist(7),'color','m');
%         title(ax,sprintf('gs %.2f',gridscore2(aco,-1)),'color','m','fontsize',18);
%         %AC CC SMOOTHED
%         ax = subplot(rows,cols,8,'Parent',m.parent); axi(end+1) = ax;
%         %ac = imgaussfilt(aco, 2,'FilterDomain','spatial'); %imshow same
%         imagesc(ax,ac);hold(ax,'on');





%imagesc(ax,imgaussfilt(rm, 1)); %imshow
% CORRELATION
%imagesc(corrcoef(rm));
%imagesc(xcorr2(rm,rm));
%a(a > 0.1) = 0.1;
%imagesc(a);



