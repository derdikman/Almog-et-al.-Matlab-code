function TabTime(fig, tab, m)
    
    gui.Window = fig;
    gui.m = m;
    gui.m.tab = tab;
    gui.m.sesh = 'before';
    gui.m.version = 'n';
    gui.m.bindeg = 360;
    gui.m.lag = 2.4;
    gui.m.binspike = 0.06;
    gui.m.overlap = 0.002;
    gui.m.sigma = 0;
    gui.m.movmean = 0;
    gui.m.bandf = -1;
    gui.m.bandw = 0;
    gui.Window.UserData.gidsFns{end+1} = @updateGids;
    
    gui.Window.UserData.nza = 34;
    
    %top box
    TimeTabBox = uix.HBoxFlex( 'Parent', tab, 'Spacing', 3);
    gui.paramPanel = uix.Panel('Parent', TimeTabBox,'Title', 'Parameters' );
    gui.viewPanel = uix.Panel('Parent', TimeTabBox,'Title', 'View','fontweight','bold',...
        'FontSize', 11, 'TitlePosition','centerTop','ForegroundColor','r');
    set( TimeTabBox, 'Widths', [200,-1] );% Adjust the main layout
    %% VIEW PANEL
    viewLayout = uix.VBox( 'Parent', gui.viewPanel,'Padding', 0, 'Spacing', 0);
    t = [];
    gui.viewContainer = uicontainer('Parent', viewLayout); t = [t -1]; %'backgroundColor','w'
    set(viewLayout, 'Heights', t);
    %% PARAMETERS PANEL
    pad = 1; spc = 1;
    paramsLayout = uix.VBox( 'Parent', gui.paramPanel,'Padding', pad,'Spacing',spc);
    %DEG
    degBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', pad, 'Spacing', spc );
    uicontrol('Parent', degBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Direction Bins(deg):');
    degBoxBoxH = uix.HBox( 'Parent', degBoxV,'Padding', pad, 'Spacing', spc );
    step = 3;step = step/(360-3); %max-min
    gui.degSlider = uicontrol( 'Style', 'slider', 'Min', 3,'Max', 360,'Parent', degBoxBoxH, ...
        'Value', 360,'SliderStep',[step 3*step],'Callback', @onDegSlide);
    gui.degEdit = uicontrol( 'Parent', degBoxBoxH, 'Style', 'edit', 'String', 360,...
        'Callback', @onDegEdit);
    addlistener(gui.degSlider,'ContinuousValueChange',@(hObject, event) onDegSliding(hObject, event));
    set( degBoxBoxH, 'Widths', [-4 -1] );
    %LAG
    lagBoxV = uix.VBox( 'Parent', paramsLayout,'Padding',pad,'Spacing',spc, 'Visible', 'on');
    uicontrol('Parent', lagBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Correlation Lag(s):');
    lagBoxH = uix.HBox( 'Parent', lagBoxV,'Padding',spc,'Spacing',pad);
    step = 0.1; step = step/(5-0.1); %max-min
    gui.lagSlider = uicontrol( 'Style', 'slider', 'Min', 0.1,'Max', 5.0, ...
        'Parent', lagBoxH,'Value', 2.4,'SliderStep',[step 3*step],'Callback', @onLagSlide);
    gui.lagEdit = uicontrol( 'Parent', lagBoxH, 'Style', 'edit', 'String', 2.4,...
        'Callback', @onLagEdit);
    addlistener(gui.lagSlider,'ContinuousValueChange',@(hObject, event) onLagSliding(hObject, event));
    set( lagBoxH, 'Widths', [-4 -1] );
    %SPIKE
    spikeBinBoxV = uix.VBox( 'Parent', paramsLayout,'Padding',pad,'Spacing',spc, 'Visible', 'on');
    uicontrol('Parent', spikeBinBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Spike Time Bin(s):');
    spikeBinBoxH = uix.HBox( 'Parent', spikeBinBoxV, ...
        'Padding',spc,'Spacing',pad);
    mnx = [.02 .2]; step = 0.02; step = step/(mnx(2)-mnx(1));%step = 0.02; step = step/(0.2-0.02); %max-min
    gui.spikeBinSlider = uicontrol( 'Style', 'slider','Min', mnx(1),'Max', mnx(2), 'SliderStep',[step 2*step],...
        'Parent', spikeBinBoxH,'Value', gui.m.binspike,'Callback', @onSpikeBinSlide);
    gui.spikeBinEdit = uicontrol( 'Parent', spikeBinBoxH, 'Style', 'edit', 'String', gui.m.binspike, ...
        'Callback', @onSpikeBinEdit);
    addlistener(gui.spikeBinSlider,'ContinuousValueChange',@(hObject, event) onSpikeSliding(hObject, event));
    set( spikeBinBoxH, 'Widths', [-4 -1] );
    %OVERLAP
        overlapBoxV = uix.VBox( 'Parent', paramsLayout,'Padding',pad,'Spacing',spc, 'Visible', 'on' );
    uicontrol('Parent', overlapBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Remove overlapping spikes [0-0.5s]:');
    overlapBoxH = uix.HBox( 'Parent', overlapBoxV,'Padding',spc,'Spacing',pad);
    mnx = [0 0.5]; step = 0.001; step = step/(mnx(2)-mnx(1)); %max-min
    gui.overlapSlider = uicontrol( 'Style', 'slider','Min', mnx(1),'Max', mnx(2), 'SliderStep',[step 2*step],...
        'Parent', overlapBoxH,'Value', gui.m.overlap,'Callback', @onOverlapSlide);
    gui.overlapEdit = uicontrol( 'Parent', overlapBoxH, 'Style', 'edit', 'String', gui.m.overlap, ...
        'Callback', @onOverlapEdit);
    addlistener(gui.overlapSlider,'ContinuousValueChange',@(hObject, event) onOverlapSliding(hObject, event));
    set( overlapBoxH, 'Widths', [-4 -1] );
     %SIGMA
    sigmaBoxV = uix.VBox( 'Parent', paramsLayout,'Padding',pad,'Spacing',spc, 'Visible', 'on' );
    uicontrol('Parent', sigmaBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Low pass train (hz)');
    sigmaBoxH = uix.HBox( 'Parent', sigmaBoxV,'Padding',spc,'Spacing',pad);
    mnx = [0 500]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
    gui.sigmaSlider = uicontrol( 'Style', 'slider','Min', mnx(1),'Max', mnx(2), 'SliderStep',[step 2*step],...
        'Parent', sigmaBoxH,'Value', gui.m.sigma,'Callback', @onSigmaSlide);
    gui.sigmaEdit = uicontrol( 'Parent', sigmaBoxH, 'Style', 'edit', 'String', gui.m.sigma, ...
        'Callback', @onSigmaEdit);
    addlistener(gui.sigmaSlider,'ContinuousValueChange',@(hObject, event) onSigmaSliding(hObject, event));
    set( sigmaBoxH, 'Widths', [-4 -1] );
    %HAMWIN
        movmeanBoxV = uix.VBox( 'Parent', paramsLayout,'Padding',pad,'Spacing',spc, 'Visible', 'on' );
    uicontrol('Parent', movmeanBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Moving avg (ms)');
    movmeanBoxH = uix.HBox( 'Parent', movmeanBoxV,'Padding',spc,'Spacing',pad);
    mnx = [0 1000]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
    gui.movmeanSlider = uicontrol( 'Style', 'slider','Min', mnx(1),'Max', mnx(2), 'SliderStep',[step 2*step],...
        'Parent', movmeanBoxH,'Value', gui.m.movmean,'Callback', @onMovmeanSlide);
    gui.movmeanEdit = uicontrol( 'Parent', movmeanBoxH, 'Style', 'edit', 'String', gui.m.movmean, ...
        'Callback', @onMovmeanEdit);
    addlistener(gui.movmeanSlider,'ContinuousValueChange',@(hObject, event) onMovmeanSliding(hObject, event));
    set( movmeanBoxH, 'Widths', [-4 -1] );
        %BANDPASS 
        bandBoxV = uix.VBox( 'Parent', paramsLayout,'Padding',pad,'Spacing',spc, 'Visible', 'on' );
    uicontrol('Parent', bandBoxV,'Style','text', 'HorizontalAlignment', 'left','String', 'Bandpass Hz: freq, bandwidth');
        %f
    bandfBoxH = uix.HBox( 'Parent', bandBoxV,'Padding',spc,'Spacing',pad);
    mnx = [-1 500]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
    gui.bandfSlider = uicontrol( 'Style', 'slider','Min', mnx(1),'Max', mnx(2), 'SliderStep',[step 2*step],...
        'Parent', bandfBoxH,'Value', gui.m.bandf,'Callback', @onBandfSlide);
    gui.bandfEdit = uicontrol( 'Parent', bandfBoxH, 'Style', 'edit', 'String', gui.m.bandf,'Callback', @onBandfEdit);
    addlistener(gui.bandfSlider,'ContinuousValueChange',@(hObject, event) onBandfSliding(hObject, event));
    set( bandfBoxH, 'Widths', [-4 -1] );
        %w
    bandwBoxH = uix.HBox( 'Parent', bandBoxV,'Padding',spc,'Spacing',pad);
    mnx = [0 300]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min
    gui.bandwSlider = uicontrol( 'Style', 'slider','Min', mnx(1),'Max', mnx(2), 'SliderStep',[step 2*step],...
        'Parent', bandwBoxH,'Value', gui.m.bandw,'Callback', @onBandwSlide);
    gui.bandwEdit = uicontrol( 'Parent', bandwBoxH, 'Style', 'edit', 'String', gui.m.bandw,'Callback', @onBandwEdit);
    addlistener(gui.bandwSlider,'ContinuousValueChange',@(hObject, event) onBandwSliding(hObject, event));
    set( bandwBoxH, 'Widths', [-4 -1] );
    
    %GROUP
    groupBoxV = uix.VBox( 'Parent', paramsLayout,'Padding',spc,'Spacing',pad);
    uicontrol('Parent', groupBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Group:');
    gui.groupListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','Parent', groupBoxV, ...
        'String', gui.Window.UserData.gids,'Value', 1,'Callback', @onGroupListSelection);
    %MID
    midBoxV = uix.VBox( 'Parent', paramsLayout,'Padding',spc,'Spacing',pad);
    uicontrol('Parent', midBoxV,'Style','text','HorizontalAlignment', ...
        'left', 'String','Muscimol Bin:');
    gui.midListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w', ...
        'Parent', midBoxV,'String', {'before';'midall';'after'}, ...
        'Value', 1,'Callback', @onMidListSelection, 'Enable', 'on');
    if ~gui.m.g(1).after.exists
        set(gui.midListBox,'String', {'before';'midall'});
    end
    %RUN
    runBoxV = uix.VBox( 'Parent', paramsLayout,'Padding',spc,'Spacing',pad);
    gui.RunButton = uicontrol( 'Style', 'PushButton','Parent', runBoxV,'String', 'Analysis!', ...
        'Callback', @onRunButton );
    gui.status = uicontrol('Parent', runBoxV,'Style','text', 'Tag', 'status',...
        'HorizontalAlignment', 'left', 'String', 'ok', 'FontSize',9);
    %% PARAMETERS PANEL TIGHTEN UP
    %Tighten UP
    t = []; htlab = 15; a = htlab + 20; 
    set(degBoxV, 'Heights', [htlab -1]); t = [t a]; % Make the list fill the space
    set(lagBoxV, 'Heights', [htlab -1]); t = [t a];% Make the list fill the space
    set(spikeBinBoxV, 'Heights', [htlab -1]); t = [t a]; % Make the list fill the space
    set(sigmaBoxV, 'Heights', [htlab -1]); t = [t a]; % Make the list fill the space
    set(overlapBoxV, 'Heights', [htlab -1]); t = [t a]; % Make the list fill the space
    set(movmeanBoxV, 'Heights', [htlab -1]); t = [t a]; % Make the list fill the space
    set(bandBoxV, 'Heights', [htlab -1 -1]); t = [t a+20]; % Make the list fill the space
    set(groupBoxV, 'Heights', [htlab -1]); t = [t 100]; % Make the list fill the space
    set(midBoxV, 'Heights', [htlab -1]); t = [t 70]; % Make the list fill the space
    set(runBoxV, 'Heights', [30 -1]); t = [t -1]; % Make the list fill the space
    
    set(paramsLayout, 'Heights', t); % Make the lists fill the space
    
    function updateGids()
        set(gui.groupListBox, 'String', gui.Window.UserData.gids);
    end
    
    %% CALLBACKS
    %DEG
    function onDegSlide( src, ~ )
        t = round2r(get(src, 'Value'),3);
        gui.m.bindeg = t;
        set(gui.degEdit, 'String', t);
        run();
    end
    function onDegEdit( src, ~ )
        t = round2r(str2double(get(src, 'String')),3);
        gui.m.bindeg = t;
        set(gui.degSlider, 'Value', t);
    end
    function onDegSliding(hObject,event)
        t = round2r(get(hObject,'Value'),3);
        set(gui.degEdit, 'String', t);
    end
    %LAG
    function onLagSlide( src, ~ )
        t = round2r(get(src, 'Value'), 0.1);
        gui.m.lag = t;
        set(gui.lagEdit, 'String', t);
        %run();
    end
    function onLagEdit( src, ~ )
        t = round2r(str2double(get(src, 'String')),0.1);
        gui.m.lag = t;
        set(gui.lagSlider, 'Value', t);
    end
    function onLagSliding(hObject,event)
        t = round2r(get(hObject,'Value'),0.1);
        set(gui.lagEdit, 'String', t);
    end
    %SPIKE
    function onSpikeBinSlide( src, ~ )
        t = round2r(get( src, 'Value'),0.02);
        gui.m.binspike = t;
        set(gui.spikeBinEdit, 'String', t);
        %run();
    end
    function onSpikeBinEdit(src, ~ )
        t = round2r(str2double(get(src, 'String')),0.02);
        gui.m.binspike =t;
        set(gui.spikeBinSlider, 'Value', t);
    end
    function onSpikeSliding(hObject,event)
        t = round(get(hObject,'Value')/0.02)*0.02;
        set(gui.spikeBinEdit, 'String', t);
    end
    %SIGMA
    function onSigmaSlide( src, ~ )
        t = round(get( src, 'Value'));
        gui.m.sigma = t;
        set(gui.sigmaEdit, 'String', t);
        run();
    end
    function onSigmaEdit( src, ~ )
        t = round(str2double(get(src, 'String')));
        gui.m.sigma =t;
        set(gui.sigmaSlider, 'Value', t);
    end
    function onSigmaSliding(hObject,event)
        t = round2r(get(hObject,'Value'),1);
        set(gui.sigmaEdit, 'String', t);
    end
    %OVERLAP
    function onOverlapSlide( src, ~ )
        t = round(get( src, 'Value'),3);
        gui.m.overlap = t;
        set(gui.overlapEdit, 'String', t);
        %run();
    end
    function onOverlapEdit( src, ~ )
        t = round(str2double(get(src, 'String')),3);
        gui.m.overlap =t;
        set(gui.overlapSlider, 'Value', t);
    end
    function onOverlapSliding(hObject,event)
        t = round2r(get(hObject,'Value'),0.001);
        set(gui.overlapEdit, 'String', t);
    end
 %HAMWIN
    function onMovmeanSlide( src, ~ )
        t = round(get( src, 'Value'));
        gui.m.movmean = t;
        set(gui.movmeanEdit, 'String', t);
        run();
    end
    function onMovmeanEdit( src, ~ )
        t = round(str2double(get(src, 'String')),3);
        gui.m.movmean =t;
        set(gui.movmeanSlider, 'Value', t);
    end
    function onMovmeanSliding(hObject,event)
        t = round2r(get(hObject,'Value'),1);
        set(gui.movmeanEdit, 'String', t);
    end
 %Band F
    function onBandfSlide( src, ~ )
        t = round(get( src, 'Value'));
        gui.m.bandf = t;
        set(gui.bandfEdit, 'String', t);
        %run();
    end
    function onBandfEdit( src, ~ )
        t = round(str2double(get(src, 'String')),3);
        gui.m.bandf =t;
        set(gui.bandfSlider, 'Value', t);
    end
    function onBandfSliding(hObject,event)
        t = round2r(get(hObject,'Value'),1);
        set(gui.bandfEdit, 'String', t);
    end
 %BAND W
    function onBandwSlide( src, ~ )
        t = round(get( src, 'Value'));
        gui.m.bandw = t;
        set(gui.bandwEdit, 'String', t);
        %run();
    end
    function onBandwEdit( src, ~ )
        t = round(str2double(get(src, 'String')),3);
        gui.m.bandw =t;
        set(gui.bandwSlider, 'Value', t);
    end
    function onBandwSliding(hObject,event)
        t = round2r(get(hObject,'Value'),1);
        set(gui.bandwEdit, 'String', t);
    end
    % GROUP FNS
    function onGroupListSelection( src, t )
        set(gui.status, 'ForegroundColor',[0,0,1]);
           t = cellstr(get( src, 'String' ));
        gui.m.gid = str2double(t(get( src, 'Value' )));
        gui.m.g = gui.m.groups{gui.m.gid};
        l = length(gui.m.g);
        set(gui.status, 'String', sprintf('group size: %d',l));
        set(gui.midListBox,'String',{'before';'midall';'after'},'Value', 1);
        %gui.m.mid = 1;
        %drawnow(); %TRY??
        run();
    end
    % MID
    function onMidListSelection( src, e )
        %gui.m.mid = get( src, 'Value' )-1;
        t = get( src, 'String');
        gui.m.sesh = str2double(t{get( src, 'Value' )});
        if isnan(gui.m.sesh )
            gui.m.sesh = t{get( src, 'Value' )};
        end
    end
    
    function onRunButton( ~, ~ )
        run();
    end % onDemoHelp
    function onExit( ~, ~ )
        % User wants to quit out of the application
        delete( gui.Window );
    end % onExit
    %% RUN
    function run()
        tic
        set(gui.status, 'ForegroundColor',[0,0,1]);
        disp(['time run g', num2str(gui.m.gid)]);
        delete(findobj(gui.m.tab,'type','axes'));
        drawnow();
        gui.m.grid_thresh = gui.Window.UserData.gridThresh;
        params = gui.m;
        set(gui.status, 'String', 'computing....');
        %enable(handles, false)
        params.sesh = params.sesh;
        params.parent = gui.viewContainer;
        params.fig = gui.Window;
        params.good = gui.m.g(gui.Window.UserData.cids{gui.m.gid});
        [~,v] = plotByDirectionMainTimeTab(params);
        set(gui.status, 'String', sprintf('%s',v));
        set(gui.status, 'ForegroundColor',[0,0,1]);
        set(gui.status, 'String', v);
        set(gui.viewPanel,'Title', v);
        'run'
        toc
        %enable(handles, true);
    end
end


