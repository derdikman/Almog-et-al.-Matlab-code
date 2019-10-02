function TabMotion(fig, tab, m)
    gui.Window = fig;
    gui.m = m;
    gui.m.speed = 200;
    gui.m.sesh = 'before';
    gui.m.cells = ones(length(gui.m.g),1);
    gui.m.pause = false;
    gui.m.stop = false;
    gui.m.delay = 0;
    gui.m.go = 1;
    gui.m.t0 = 1;
    gui.m.rmmode = false;
    gui.m.lastChecked = 1;
    gui.m.timestep = 1;
    gui.m.win = 500;
    gui.m.showTrain = false;
    gui.m.drift = 60;
    gui.m.showDrift = false;
    %gui.m.movie = VideoWriter('gridmotion.avi'); open(gui.m.movie);
    gui.m.record = false;
    updateM();
    maxx  = 0; maxy = 0;
    
    %TOP UI
    AniTabBox = uix.HBoxFlex( 'Parent', tab, 'Spacing', 3);
    controlPanel = uix.Panel('Parent', AniTabBox,'Title', '' );
    gui.viewPanel = uix.Panel('Parent', AniTabBox,'Title', 'Ready','fontweight','bold',...
        'FontSize', 11, 'TitlePosition','centerTop','ForegroundColor','b');
    set( AniTabBox, 'Widths', [200 -1] );
    
    %% VIEW PANEL
    viewBoxV = uix.VBoxFlex( 'Parent', gui.viewPanel, 'Spacing', 3,'backgroundcolor','k');
    gui.viewContainer  = uicontainer('Parent', viewBoxV,'backgroundcolor','k'); %gui.viewPanel
    gui.viewContainer2 = uicontainer('Parent', viewBoxV,'backgroundcolor','k');
    set( viewBoxV, 'Heights', [-7 -1] );
    %% CONTROL PANEL
    controlBoxV = uix.VBox( 'Parent', controlPanel, 'Spacing', 1);
    t = [];
    pad = 1; spc = 1; sld = 35;
    %% GROUP
    groupBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', pad, 'Spacing', spc ); t = [t -1];
    uicontrol('Style','text','Parent', groupBoxV,'HorizontalAlignment', 'left','String', 'Group:');
    gui.groupListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','Parent', groupBoxV, ...
        'String', num2str((1:length(gui.m.groups))'),'Value', 1,'Callback', @onGroupListSelection);
    set( groupBoxV, 'Heights', [20 -1] );
    %% MID
    
    midBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', pad, 'Spacing', spc ); t = [t 70];
    uicontrol('Parent', midBoxV,'Style','text','HorizontalAlignment', ...
        'left', 'String','Muscimol Bin:');
    gui.midListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','String', {'before';'midall';'after'}, ...
        'Parent', midBoxV, 'Value', 1,'Callback', @onMidListSelection, 'Enable', 'on');
    if ~gui.m.g(1).after.exists
        set(gui.midListBox,'String', {'before';'midall'});
    end
    set( midBoxV, 'Heights', [20 -1] );
    %% CHECK BOXES
    gui.checkBoxV = uix.VBox('Parent', controlBoxV); t = [t 200];
    uicontrol('Style','text','Parent', gui.checkBoxV, 'String', 'Show Cells:',...
        'HorizontalAlignment', 'left');
    gui.cbh = zeros(length(gui.m.g),1);
    for i = 1:length(gui.m.g)
        c = getSesh(i);
        s = sprintf('i%d tet%d grid%.1f', gui.m.g(i).ind, gui.m.g(i).tet, c.gridscore);
        gui.cbh(i) = uicontrol('Style','checkbox','String', s,'Value',1,...
                'Parent', gui.checkBoxV, 'Callback',{@checkBoxCallback,i},...
                'ForegroundColor',gui.m.colors(i,:), 'FontSize', 7, 'fontweight', 'bold');
    end
    %% DRIFT
    driftBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 0, 'Spacing',0 , 'Visible', 'on'); t = [t sld];
    driftBoxH1 = uix.HBox( 'Parent', driftBoxV,'Padding', 0, 'Spacing', 0 );
    gui.showDriftCh = uicontrol('Style','checkbox','Value',gui.m.showTrain,'String', 'Drift Window(s):',...
          'Parent', driftBoxH1, 'Callback',@OnShowDriftCh); %'FontSize', 10, 'fontweight', 'bold'
    %uicontrol('Parent', driftBoxH1,'Style','text','String', 'Train Window(s):', 'HorizontalAlignment', 'left');
    driftBoxH = uix.HBox( 'Parent', driftBoxV,'Padding', 0, 'Spacing', 0 );
    mnx = [10 600]; step = 10; step = step/(mnx(2)-mnx(1)); %max-min 
    gui.driftSlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
        'Parent', driftBoxH,'Value', gui.m.drift,'SliderStep',[step 3*step],'Callback', @onDriftSlide);        
    gui.driftEdit = uicontrol( 'Parent', driftBoxH, 'Style', 'edit', 'String', gui.m.drift,...
        'Callback', @onDriftEdit);
    addlistener(gui.driftSlider,'ContinuousValueChange',@(hObject, event) onDriftSliding(hObject, event));
    set( driftBoxH, 'Widths', [-4 -1] );
    %% GO
    goBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 0, 'Spacing', 0 , 'Visible', 'on'); t = [t sld];      
    gui.goLabel = uicontrol('Parent', goBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', sprintf('Go to (min %d): ',floor(gui.m.c.pt(end)/60)));
    goBoxH = uix.HBox( 'Parent', goBoxV,'Padding', 0, 'Spacing', 0 );
    mnx = [0 floor(double(gui.m.c.pt(end))/60)]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min 
    gui.goSlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
        'Parent', goBoxH,'Value', 0,'SliderStep',[step 2*step],'Callback', @onGoSlide);        
    gui.goEdit = uicontrol( 'Parent', goBoxH, 'Style', 'edit', 'String', 0,...
        'Callback', @onGoEdit);
    addlistener(gui.goSlider,'ContinuousValueChange',@(hObject, event) onGoSliding(hObject, event));
    set( goBoxH, 'Widths', [-4 -1] );
    %% T0
    t0BoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 0, 'Spacing', 0 , 'Visible', 'on'); t = [t sld];      
    gui.t0Label = uicontrol('Parent', t0BoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', sprintf('Show from t= (min <%d): ',floor(gui.m.c.pt(end)/60)));
    t0BoxH = uix.HBox( 'Parent', t0BoxV,'Padding', 0, 'Spacing',0 );
    mnx = [0 floor(double(gui.m.c.pt(end))/60)]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min 
    gui.t0Slider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
        'Parent', t0BoxH,'Value', 0,'SliderStep',[step 2*step],'Callback', @onT0Slide);        
    gui.t0Edit = uicontrol( 'Parent', t0BoxH, 'Style', 'edit', 'String', 0,...
        'Callback', @onT0Edit);
    addlistener(gui.t0Slider,'ContinuousValueChange',@(hObject, event) onT0Sliding(hObject, event));
    set( t0BoxH, 'Widths', [-4 -1] );   
    %% TIME
    timeBoxH = uix.VBox( 'Parent', controlBoxV,'Padding', 0, 'Spacing', 0 , 'Visible', 'on'); t = [t sld];      
    uicontrol('Parent', timeBoxH,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Time Step:');
    timeBoxH = uix.HBox( 'Parent', timeBoxH,'Padding', 0, 'Spacing', 0 );
    mnx = [1 1000]; step = 1; step = step/(mnx(2)-mnx(1)); %max-min 
    gui.timeSlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
        'Parent', timeBoxH,'Value', gui.m.timestep,'SliderStep',[step 10*step],'Callback', @onTimeSlide);        
    gui.timeEdit = uicontrol( 'Parent', timeBoxH, 'Style', 'edit', 'String', gui.m.timestep,...
        'Callback', @onTimeEdit);
    addlistener(gui.timeSlider,'ContinuousValueChange',@(hObject, event) onTimeSliding(hObject, event));
    set( timeBoxH, 'Widths', [-4 -1] );
    %% SPEED
    speedBoxH = uix.VBox( 'Parent', controlBoxV,'Padding', 0, 'Spacing', 0 , 'Visible', 'on'); t = [t sld];      
    uicontrol('Parent', speedBoxH,'Style','text','String', 'Speed (lightyears):',...
        'HorizontalAlignment', 'left');
    speedBoxH = uix.HBox( 'Parent', speedBoxH,'Padding', 0, 'Spacing', 0 );
    mnx = [1 2000]; step = 10; step = step/(mnx(2)-mnx(1)); %max-min 
    gui.speedSlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
        'Parent', speedBoxH,'Value', gui.m.speed,'SliderStep',[step 3*step],'Callback', @onSpeedSlide);        
    gui.speedEdit = uicontrol( 'Parent', speedBoxH, 'Style', 'edit', 'String', gui.m.speed,...
        'Callback', @onSpeedEdit);
    addlistener(gui.speedSlider,'ContinuousValueChange',@(hObject, event) onSpeedSliding(hObject, event));
    set( speedBoxH, 'Widths', [-4 -1] );
     %% DELAY
    delayBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 0, 'Spacing', 0 , 'Visible', 'on'); t = [t sld];      
    gui.delayLabel = uicontrol('Parent', delayBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', sprintf('Delay (s):%s',''));
    delayBoxH = uix.HBox( 'Parent', delayBoxV,'Padding', 0, 'Spacing', 0 );
    mnx = [-3 3]; step = 0.02; step = step/(mnx(2)-mnx(1)); %max-min 
    gui.delaySlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
        'Parent', delayBoxH,'Value', gui.m.delay,'SliderStep',[step 5*step],'Callback', @onDelaySlide);        
    gui.delayEdit = uicontrol( 'Parent', delayBoxH, 'Style', 'edit', 'String', gui.m.delay,...
        'Callback', @onDelayEdit);
    addlistener(gui.delaySlider,'ContinuousValueChange',@(hObject, event) onDelaySliding(hObject, event));
    set( delayBoxH, 'Widths', [-4 -1] );   
    %% TRAIN DISPLAY WIN
    winBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 0, 'Spacing', 0 , 'Visible', 'on'); t = [t sld];
    winBoxH1 = uix.HBox( 'Parent', winBoxV,'Padding', 0, 'Spacing', 0 );
    gui.showTrainCh = uicontrol('Style','checkbox','Value',gui.m.showTrain,'String', 'Train Length(s):',...
          'Parent', winBoxH1, 'Callback',@OnShowTrainCh); %'FontSize', 10, 'fontweight', 'bold'
    %uicontrol('Parent', winBoxH1,'Style','text','String', 'Train Window(s):', 'HorizontalAlignment', 'left');
    winBoxH = uix.HBox( 'Parent', winBoxV,'Padding', 0, 'Spacing', 0 );
    mnx = [50 20000]; step = 50; step = step/(mnx(2)-mnx(1)); %max-min 
    gui.winSlider = uicontrol( 'Style', 'slider', 'Min', mnx(1),'Max', mnx(2), ...
        'Parent', winBoxH,'Value', gui.m.win,'SliderStep',[step 3*step],'Callback', @onWinSlide);        
    gui.winEdit = uicontrol( 'Parent', winBoxH, 'Style', 'edit', 'String', gui.m.win,...
        'Callback', @onWinEdit);
    addlistener(gui.winSlider,'ContinuousValueChange',@(hObject, event) onWinSliding(hObject, event));
    set( winBoxH, 'Widths', [-4 -1] );
    %% BUTTONS
    buttonBoxH = uix.HBox( 'Parent', controlBoxV, 'Spacing', 3); t = [t 50];
    gui.RunButton = uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'Start','Callback', @run );
    gui.PauseButton = uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'Pause','Callback', @onPause,'Enable','On');
    gui.StopButton = uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'Stop','Callback', @onStop,'Enable','On');
    gui.RecordButton = uicontrol( 'Style', 'PushButton','Parent', buttonBoxH,...
        'String', 'Record','Callback', @onRecord,'Enable','On','backgroundColor','w');
    %STATUS
    %gui.status = uicontrol('Parent', controlBoxV,'Style','text', 'Tag', 'status','String', 'listo',...
    %'HorizontalAlignment', 'left', 'FontSize',9,'foregroundcolor','b'); t = [t 30];
    %% FINAL HEIGHTS
    set(controlBoxV, 'Heights', t );
    
    %% CALLBACKS
    % CHECKBOXES
    function checkBoxCallback(hObject,eventData,id)
        gui.m.cells(id) = get(hObject,'Value'); %[id get(hObject,'Value')]
        %if get(hObject,'Value')
        gui.m.lastChecked = id;
        %end
        setPTitle();
    end
    % GROUP
    function onGroupListSelection( src, ~ )
        %set(gui.status, 'ForegroundColor',[0,0,1]);
        onStop(-1,-1);
        set(gui.viewPanel,'ForegroundColor',[0,0,1], 'Title', '');
        t = cellstr(get( src, 'String' ));
        gui.m.gid = str2double(t(get( src, 'Value' )));
        gui.m.g = gui.m.groups{gui.m.gid};
        gui.m.sesh = 'before';
        gui.m.lastChecked = 1;
        updateM();
        disp(['G: ', num2str(gui.m.gid),' ', num2str(length(gui.m.g))]);
        %set(gui.viewPanel,'Title', sprintf('group size: %d',length(gui.m.g)));
        %set(gui.status, 'String', sprintf('group size: %d',l));cellstr((num2str(1:length(gui.m.g(1).middle),'%d'))');
        set(gui.midListBox,'String',{'before';'midall';'after'},'Value', 1);
        if ~gui.m.g(1).after.exists
            set(gui.midListBox,'String', {'before';'midall'});
        end
        %reset checkboxes
        delete(findobj(gui.checkBoxV.Children,'Style','checkbox'));
        gui.cbh = zeros(length(gui.m.g),1);
        gui.m.cells = ones(length(gui.m.g),1);
        for i = 1:length(gui.m.g)
            c = getSesh(i);
            s = sprintf('i%d tet%d grid%.1f', gui.m.g(i).ind, gui.m.g(i).tet, c.gridscore);
            gui.cbh(i) = uicontrol('Style','checkbox','String', s,'Value',1,...
                'Parent', gui.checkBoxV, 'Callback',{@checkBoxCallback,i},...
                'ForegroundColor',gui.m.colors(i,:), 'FontSize', 9, 'fontweight', 'bold');
        end
        %[~,s] = getSesh(1);
        %set(gui.viewPanel,'Title', s);
        setPTitle();
        gui.m.go = 1; %$$$$$$$$$$$$$$$
        gui.m.t0 = 1;
        maxx = max(gui.m.c.px);
        maxy = max(gui.m.c.py);
        %run();  %REDO THIS
    end
    % MID 
    function onMidListSelection( src, ~ )
        %gui.m.mid = get( src, 'Value' )-1;
        gui.m.stop = true;
        t = get( src, 'String');
        gui.m.sesh = str2double(t{get( src, 'Value' )});
        if isnan(gui.m.sesh)
            gui.m.sesh = t{get(src,'Value')};
        end
        %disp('1 MID SELCTION');
        updateM();
        %checkboxes
        for i = 1:length(gui.m.g)
            c = getSesh(i);
            s = sprintf('i%d tet%d grid%.1f', gui.m.g(i).ind, gui.m.g(i).tet, c.gridscore);
            set(gui.cbh(i),'String', s);
        end
        setPTitle();
        %run(); %REDO THIS
        %}
    end
    % STOP
    function onStop(~, ~)
        %disp('STOP');
        gui.m.stop = true;
        gui.m.go = 1;
        gui.m.t0 = 1;
    end
    % PAUSE
    function onPause(~, ~)
        %delete(findobj(gui.Window.Children.Children(gui.Window.Children.Selection),...
        %'type','axes'));
        if gui.m.pause
            set(gui.PauseButton, 'String', 'Pause')
        else
            set(gui.PauseButton, 'String', 'Resume')
        end
        drawnow
        gui.m.pause = ~gui.m.pause;
    end
    % Record
    function onRecord(~, ~)
        if gui.m.record %stop
            %figure;
            %movie(gui.m.movie);
            close(gui.m.movie);
            set(gui.RecordButton, 'String', 'Record','backgroundColor','w')
        else %start
            %gui.m.movie = VideoWriter(sprintf('Grid_Motion_%s.mp4',datestr(now,'yyyy.mm.dd_HH.MM.SS'))); %','
            gui.m.movie = VideoWriter(sprintf('Grid_Motion_%s.avi',datestr(now,'yyyy.mm.dd_HH.MM.SS'))); 
            %gui.m.movie.Quality = 100;
            open(gui.m.movie);
            set(gui.RecordButton, 'String', 'Record','backgroundColor','r');
        end
        gui.m.record = ~gui.m.record ;
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
    % Go
    function onGoSlide( src, ~ )
        c = getSesh(1);
        t = round2r(get(src, 'Value'), 1);
        gui.m.go = round(t*60/.02);
        set(gui.goEdit, 'String', t);
    end
    function onGoEdit(src,~)
        c = getSesh(1);
        t = round2r(str2double(get(src,'String')), 1); %str2double(get(src,'String'))
        gui.m.go = round(t*60/.02);
        set(gui.goSlider, 'Value', t);
    end 
    function onGoSliding(hObject,~)
        t = round2r(get(hObject,'Value'),1);
        set(gui.goEdit, 'String', t);
    end
    % T0
    function onT0Slide( src, ~ )
        c = getSesh(1);
        t = round2r(get(src, 'Value'), 1);
        gui.m.t0 = round(t*60/.02)+1;
        set(gui.t0Edit, 'String', t);
    end
    function onT0Edit(src,~)
        c = getSesh(1);
        t = round2r(str2double(get(src,'String')), 1); %str2double(get(src,'String'))
        gui.m.t0 = round(t*60/.02)+1;
        set(gui.t0Slider, 'Value', t);
    end 
    function onT0Sliding(hObject,~)
        t = round2r(get(hObject,'Value'),1);
        set(gui.t0Edit, 'String', t);
    end     
    % TIMESTEP
    function onTimeSlide( src, ~ )
        t = round2r(get(src, 'Value'), 1);
        gui.m.timestep = t;
        set(gui.timeEdit, 'String', t);
        %run();
    end
    function onTimeEdit( src, ~ )
        t = round2r(str2double(get(src, 'String')),1);
        gui.m.timestep = t;
        set(gui.timeSlider, 'Value', t);
    end
    function onTimeSliding(hObject,event)
        t = round2r(get(hObject,'Value'),1);
        set(gui.timeEdit, 'String', t);
    end 
    % SPEED
    function onSpeedSlide( src, ~ )
        t = round2r(get(src, 'Value'), 1);
        gui.m.speed = t;
        set(gui.speedEdit, 'String', t);
        %run();
    end
    function onSpeedEdit( src, ~ )
        t = round2r(str2double(get(src, 'String')),1);
        gui.m.speed = t;
        set(gui.speedSlider, 'Value', t);
    end
    function onSpeedSliding(hObject,event)
        t = round2r(get(hObject,'Value'),1);
        set(gui.speedEdit, 'String', t);
    end
    % WINDOW
    function OnShowTrainCh(src, ~ )
        gui.m.showTrain = get(src, 'Value');
    end
    function onWinSlide( src, ~ )
        t = round2r(get(src, 'Value'), 100);
        gui.m.win = t;
        set(gui.winEdit, 'String', t);
        gui.m.lx = linspace(0,100,gui.m.win);
        gui.m.ly = repmat((length(gui.m.g):-1:1)',1,gui.m.win)*2 + max(gui.m.c.py) + 2;
        %run();
    end
    function onWinEdit( src, ~ )
        t = round2r(str2double(get(src, 'String')),100);
        gui.m.win = t;
        gui.m.lx = linspace(0,100,gui.m.win);
        gui.m.ly = repmat((length(gui.m.g):-1:1)',1,gui.m.win)*2 + max(gui.m.c.py) + 2;
        set(gui.winSlider, 'Value', t);
    end
    function onWinSliding(hObject,event)
        t = round2r(get(hObject,'Value'),100);
        set(gui.winEdit, 'String', t);
    end
    % DRIFT
    function OnShowDriftCh(src, ~ )
        gui.m.showDrift = get(src, 'Value');
        gui.m.rmmode = get(src, 'Value');
        setPTitle();
    end
    function onDriftSlide( src, ~ )
        t = round2r(get(src, 'Value'), 1);
        gui.m.drift = t;
        set(gui.driftEdit, 'String', t);
        setPTitle();
    end
    function onDriftEdit( src, ~ )
        t = round2r(str2double(get(src, 'String')),1);
        gui.m.drift = t;
        set(gui.driftSlider, 'Value', t);
        setPTitle();
    end
    function onDriftSliding(hObject,event)
        t = round2r(get(hObject,'Value'),1);
        set(gui.driftEdit, 'String', t);
    end 
    %% UPDATE M
    function updateM()
        %disp('2 UPDATE M');
        gui.m.g = gui.m.groups{gui.m.gid};
        gui.m.pause = false;
        gui.m.stop = false; 
        gui.m.c = getSesh(gui.m.lastChecked);
        if gui.m.go > length(gui.m.c.pt)
          gui.m.go = max(gui.m.c.pt(1),length(gui.m.c.pt) - round(1/0.02)); %set it back 1 min from end
        end
        if gui.m.t0 > gui.m.go
          gui.m.t0 = gui.m.go;
        end
        if isfield(gui, 'goSlider')
            mnx = [0 floor((gui.m.c.pt(end)-gui.m.c.pt(1))/60)]; mnx = double(mnx);
            step = 1; step = step/(mnx(2)-mnx(1)); %max-min 
            set(gui.goSlider,'Min', mnx(1),'Max', mnx(2),'Value', 0,'SliderStep',[step 2*step]);
            set(gui.goEdit,'String',0);
            set(gui.goLabel,'String', sprintf('Go to (min %d): ',mnx(2)));
        end
        makeColors();
        train();
        gui.m.lx = linspace(0,100,gui.m.win);
        gui.m.ly = repmat((length(gui.m.g):-1:1)',1,gui.m.win)*2 + max(gui.m.c.py) + 2;
    end
    %% RUN %%
    function run(~,~)
        %disp('3 RUN');
        onStop(-1,-1);
        pause(0.2);
        gui.m.stop = false;
        delete(findobj(tab,'type','axes'));
        [~,s] = getSesh(1);
        set(gui.viewPanel,'Title', s);
        gui.m.parent = gui.viewContainer;
        gui.m.go = 1;%gui.m.go + 1; %so not 0;
        %train();
        z = -1;
        runPrivate(gui.m, z);
        %gui.m.rmmode = true;
    end %end run()
    
    function runPrivate(m, z)
        %set(gui.status, 'String', sprintf('group(%d) session(%s)',gui.m.gid, t));
        if z == -1
            z = length(gui.m.c.pt);
        end
        %DRAW LOOP
        i = 1;
        %tic
        %START WHILE LOOP
        c = gui.m.c;
        mf = 1;
        xoff = 0; yoff = 0;%xoff = 5; yoff = 5;
        c.px = c.px + xoff;c.py = c.py + yoff;
        maxx = max(c.px);
        maxy = max(c.py);
        %disp('4 RUN PRIVATE');
        ax = axes('Parent',gui.m.parent);
        ax2 = axes('Parent',gui.viewContainer2);
        while gui.m.go <= length(gui.m.c.pt)
            while(gui.m.pause)
                pause(1/5);
            end
            if gui.m.go > length(c.pt)
                gui.m.go = max(c.pt(1),length(c.pt) - round(1/0.02)); %set it back 1 min from end
            end
            if gui.m.t0 > gui.m.go
                gui.m.t0 = gui.m.go;
            end
            if mod(gui.m.go,2/0.02) == 0
                set(gui.goSlider, 'Value', round(gui.m.go/(60/0.02)));
                set(gui.goEdit, 'String', num2str(round(gui.m.go/(60/0.02))));   
            end
            %xlim(ax,[0 120]); ylim(ax,[0 120]);
            i = gui.m.go;
            i1 = 1;
            delay = gui.m.delay;
            %Drift index    
            if gui.m.showDrift || gui.m.rmmode
                i1 = max(i1, gui.m.go - round(gui.m.drift/0.02));
            end
            %drift indices
            indwin = round(i1:i + delay);
            indwin = indwin(indwin>0); indwin = indwin(indwin <= length(gui.m.c.pt));
                        %RM Mode
            if gui.m.rmmode
                %if(indwin > 1000)
                pxi1 = c.px(indwin); pyi1 = c.py(indwin); pti1 = c.pt(indwin);
                %MAKE J
                sxi1 = double(pxi1(gui.m.train(gui.m.lastChecked,i1:i))); 
                syi1 = double(pyi1(gui.m.train(gui.m.lastChecked,i1:i)));
                sti1 = double(pti1(gui.m.train(gui.m.lastChecked,i1:i)));
                rm = createRateMap(pxi1, pyi1, pti1, sxi1, syi1, sti1, true, 100); 
                imagesc(ax, imgaussfilt(rm,2,'FilterDomain','spatial'));
                %title(ax, sprintf('c%d(gsc%.1f)', c1.ind,round(c1.gridscore,1)));
                axis(ax,'equal'); axis(ax,'tight'); set(ax,'ydir','normal','xticklabel',[],'yticklabel',[]);
                %xlabel(az,sprintf('Max %.0fHz',c1.max_r));
                colormap(ax,'jet');
            else
                %SPIKE PLOT MODE
                ind = round((gui.m.t0:i) + delay); 
                ind = ind(ind>0); ind = ind(ind <= length(gui.m.c.pt));
                tx = c.px(ind); ty = c.py(ind); %txi1 = c.px(indwin); tyi1 = c.py(indwin);
                %PLOTTING
                cla(ax);
                plot(ax,c.px(gui.m.t0:i),c.py(gui.m.t0:i),'Color',[0.3 0.3 0.3],'linewidth', 1);
                hold(ax, 'on');hold(ax2, 'on');
                xlim(ax, [0 maxx+20]); ylim(ax, [0  maxy+20]); xlim(ax2, [0 maxx+20]); ylim(ax2, [0  2*length(gui.m.g)+4]);
                axis(ax,'square');
                %ax.XLim
                %plot square
                %plot(ax,[xoff,xoff,maxx,maxx,xoff],[yoff,maxy+1,maxy+1,yoff,yoff],'w','linewidth', 2);
                axis(ax,'off');axis(ax2,'off');
                %l = randperm(length(gui.m.g));
                %SHOW SPIKES
                cla(ax2);
                for j = 1:length(gui.m.g)
                    %if train(j,i) && m.cells(j)
                    if gui.m.cells(j)
                        %plot(ax, double(tx(gui.m.train(j,1:i))), double(ty(gui.m.train(j,1:i)))+0.3*j,...
                        plot(ax, double(tx(gui.m.train(j,i-length(tx)+1:i))), double(ty(gui.m.train(j,i-length(tx)+1:i)))+0.3*j,...
                            '.','Color', gui.m.colors(j,:),'MarkerSize',14); %'markerfacecolor',gui.m.colors(j,:)
                        %SHOW DRIFT WINDOW
                        % plot(ax, double(txi1(gui.m.train(j,i1:i))), double(tyi1(gui.m.train(j,i1:i)))+0.3*j,...
                         %   'o','Color', gui.m.colors(j,:),'MarkerSize',10); %'markerfacecolor',gui.m.colors(j,:)
                        %plot(ax,lx(gui.m.train(j,1:i)),ly(j,gui.m.train(j,1:i)),'.','Color', colors(j,:))
                        %SHOW SPIKETRAIN ABOVE
                        wine = max(1,i-gui.m.win+1);
                        %plot(ax,lx(gui.m.train(j,wine:i)),ly(j,gui.m.train(j,wine:i)),'.','Color', colors(j,:));
                        if(gui.m.showTrain)
                            %plot(ax,gui.m.lx(gui.m.t0:i-wine+1)+xoff,gui.m.ly(j,gui.m.t0:i-wine+1)+gui.m.train(j,wine:i)*1.9+yoff,...
                                %'Color', gui.m.colors(j,:));
                            %plot(ax,gui.m.lx(gui.m.t0:i-wine+1)+xoff,gui.m.ly(j,gui.m.t0:i-wine+1)+yoff,'Color', [0.2 0.2 0.2]);
                            plot(ax2,gui.m.lx(gui.m.t0:i-wine+1)+xoff,...
                                gui.m.ly(j,gui.m.t0:i-wine+1)+gui.m.train(j,wine:i)*1.9+yoff-max(gui.m.c.py)-2,...
                                'Color', gui.m.colors(j,:),'linewidth', 3);
                            plot(ax2,gui.m.lx(gui.m.t0:i-wine+1)+xoff,gui.m.ly(j,gui.m.t0:i-wine+1)+yoff-max(gui.m.c.py)-2,...
                                'Color', [0.2 0.2 0.2],'linewidth', 3);
                        end
                    end
                    %axis(ax2,'off');
                end
            end %END RM MODE
            text(ax, double(maxx+2), double(maxy+5+yoff),sprintf('%.2f',c.pt(double(i))-c.pt(1)),...
                'fontsize',14,'color','w');
            hold(ax, 'off'); hold(ax2, 'off');%need??
            
            drawnow();
            if(gui.m.record)
                %gui.m.movie(mf) = getframe(ax);
                frame = getframe(ax);
                writeVideo(gui.m.movie,frame);
                mf = mf+1;
            end
            %delete(ax);
            pause(1/gui.m.speed);
            gui.m.go = gui.m.go + gui.m.timestep;
            if gui.m.stop
                %disp('stop and add length');
                gui.m.go = length(gui.m.c.pt) + 1;
            end
        end
        if gui.m.stop
            gui.m.stop = false;
            cla(ax);cla(ax2);
            %disp('stop and clear after loop');
        end
        gui.m.go = 1;
        %toc
    end%runPriv()
    
    function setPTitle()
        if gui.m.rmmode 
            s = sprintf('Cell %d Drift %ds',gui.m.lastChecked, gui.m.drift);
        else
            [~,s] = getSesh(1);
        end
        set(gui.viewPanel,'Title', s);
    end
    
    %% UTIL
    function makeColors()
        gui.m.colors = jet(length(gui.m.g)); %jet hsv
    end
    
    function train
        train = [];%zeros(length(m.g), length(c.pt));
        c = [];
        for i = 1:length(gui.m.g)
            c = getSesh(i);
            if strcmp(gui.m.sesh,'before')
                edges = [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf];
            elseif strcmp(gui.m.sesh,'midall')
                %edges = [-Inf, mean([c.pt(2:end); c.pt(1:end-1)]), +Inf];
                edges = [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf];
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
    
    function [c s] = getSesh(i)
        if strcmp(gui.m.sesh,'before')
            c = gui.m.g(i).before; s = 'before';
    elseif strcmp(gui.m.sesh,'midall')
            c = gui.m.g(i).midall; s = 'midall';
        elseif strcmp(gui.m.sesh,'after')
            c = gui.m.g(i).after; s = 'after';
        else
            c = gui.m.g(i).middle{gui.m.sesh}; s = ['mid' num2str(gui.m.sesh)];
        end
        s = sprintf('Group %d Session %s',gui.m.gid, s);
    end
end




function test(c, z)
    i = 0;
    figure;
    hold off
    tic
    while i < z
        i = i+1;
        plot(c.px(1:i),c.py(1:i),'k.');
        %plot(c{1}.before.px(i),c{1}.before.py(i),'k.');
        xlim([0 120]);ylim([0 120]);
        text(110,110,num2str(i),'fontsize',14,'color','k');
        text(110,110,num2str(i),'fontsize',14,'color','w');
        drawnow();
        pause(1/20000);
    end; toc
end















