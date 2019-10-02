function TabGroup(fig, tab, m)

    gui.Window = fig; 
    gui.m = m; 
    gui.m.sesh = 'before';
    goodGroups();
    % Create the panels
    %TOP UI
    DriftTabBox = uix.HBoxFlex( 'Parent', tab, 'Spacing', 3);
    gui.paramPanel = uix.Panel('Parent', DriftTabBox,'Title', 'Parameters' );
    gui.viewPanel = uix.Panel('Parent', DriftTabBox,'Title', 'View','fontweight','bold',...
        'FontSize', 11, 'TitlePosition','centerTop','ForegroundColor','r');
    set( DriftTabBox, 'Widths', [200,-1] );% Adjust the main layout
    %% VIEW PANEL
    viewLayout = uix.VBox( 'Parent', gui.viewPanel,'Padding', 0, 'Spacing', 0);
    t = [];
    gui.viewContainer = uicontainer('Parent', viewLayout); t = [t -1]; %'backgroundColor','w'
    set(viewLayout, 'Heights', t);
    %% PARAMETERS PANEL   
    paramsLayout = uix.VBox( 'Parent', gui.paramPanel,'Padding', 3, 'Spacing', 3 );
    %% GRID
    gridBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3, 'Visible', 'On');
    uicontrol('Parent', gridBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Gridscore Before > Threshold:');
    gridBoxH = uix.HBox( 'Parent', gridBoxV,'Padding', 3, 'Spacing', 3 );
    step = .1; mn = -2; mx = 2; step = step/(mx-mn); %max-min
    gui.gridSlider = uicontrol( 'Style', 'slider', 'Min', mn,'Max', mx,'Parent', gridBoxH, ...
        'Value', gui.Window.UserData.gridThresh,'SliderStep',[step 2*step],'Callback', @onGridSlide);
    gui.gridEdit = uicontrol( 'Parent', gridBoxH, 'Style', 'edit', 'String', gui.Window.UserData.gridThresh,...
        'Callback', @onGridEdit);
    addlistener(gui.gridSlider,'ContinuousValueChange',@(hObject, event) onGridSliding(hObject, event));
    set( gridBoxH, 'Widths', [-4 -1] );
    %% GRID MID
    gridMidBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3, 'Visible', 'On');
    uicontrol('Parent', gridMidBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Gridscore Middle < Threshold:');
    gridMidBoxH = uix.HBox( 'Parent', gridMidBoxV,'Padding', 3, 'Spacing', 3 );
    step = .1; mn = -2; mx = 2; step = step/(mx-mn); %max-min
    gui.gridMidSlider = uicontrol( 'Style', 'slider', 'Min', mn,'Max', mx,'Parent', gridMidBoxH, ...
        'Value', gui.Window.UserData.gridThreshMid,'SliderStep',[step 2*step],'Callback', @onGridMidSlide);
    gui.gridMidEdit = uicontrol( 'Parent', gridMidBoxH, 'Style', 'edit', 'String', gui.Window.UserData.gridThreshMid,...
        'Callback', @onGridMidEdit);
    addlistener(gui.gridMidSlider,'ContinuousValueChange',@(hObject, event) onGridMidSliding(hObject, event));
    set( gridMidBoxH, 'Widths', [-4 -1] );
    gui.Window.UserData.gridThreshMid = 0.5;
    %% LAG
    lagBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 , 'Visible', 'off');       
    uicontrol('Parent', lagBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Correlation Lag(s):');
    lagBoxH = uix.HBox( 'Parent', lagBoxV,'Padding', 3, 'Spacing', 3 );
    step = 0.1; step = step/(5-0.1); %max-min 
    gui.lagSlider = uicontrol( 'Style', 'slider', 'Min', 0.1,'Max', 5.0, ...
        'Parent', lagBoxH,'Value', 2.4,'SliderStep',[step 3*step],'Callback', @onLagSlide);        
    gui.lagEdit = uicontrol( 'Parent', lagBoxH, 'Style', 'edit', 'String', 2.4,...
        'Callback', @onLagEdit);
    addlistener(gui.lagSlider,'ContinuousValueChange',@(hObject, event) onLagSliding(hObject, event));
    set( lagBoxH, 'Widths', [-4 -1] );
    %% SPIKE
    spikeBinBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3, 'Visible', 'off');
    uicontrol('Parent', spikeBinBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Spike Time Bin(s):');
    spikeBinBoxH = uix.HBox( 'Parent', spikeBinBoxV, ...
        'Padding', 3, 'Spacing', 3 );
    step = 0.02; step = step/(0.2-0.02); %max-min
    gui.spikeBinSlider = uicontrol( 'Style', 'slider','Min', 0.02,'Max', 0.2, ...
        'Parent', spikeBinBoxH,'Value', 0.06,'SliderStep',[step 3*step],'Callback', @onSpikeBinSlide);
    gui.spikeBinEdit = uicontrol( 'Parent', spikeBinBoxH, 'Style', 'edit', 'String', 0.06, ...
        'Callback', @onSpikeBinEdit);
    addlistener(gui.spikeBinSlider,'ContinuousValueChange',@(hObject, event) onSpikeSliding(hObject, event));
    set( spikeBinBoxH, 'Widths', [-4 -1] );
    %% SIGMA
    sigmaBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3, 'Visible', 'off' );
    uicontrol('Parent', sigmaBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Gaussian(sigma):');
    sigmaBoxH = uix.HBox( 'Parent', sigmaBoxV,'Padding', 3, 'Spacing', 3 );
    step = 1; step = step/(7); %max-min 
    gui.sigmaSlider = uicontrol( 'Style', 'slider','Min', 0,'Max', 7, 'SliderStep',[step 3*step],...
        'Parent', sigmaBoxH,'Value', 1,'Callback', @onSigmaSlide);
    gui.sigmaEdit = uicontrol( 'Parent', sigmaBoxH, 'Style', 'edit', 'String', 1, ...
        'Callback', @onSigmaEdit);
    addlistener(gui.sigmaSlider,'ContinuousValueChange',@(hObject, event) onSigmaSliding(hObject, event));
    set( sigmaBoxH, 'Widths', [-4 -1] );
    %% GROUP
    groupBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    uicontrol('Parent', groupBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Group:'); %num2str((1:length(gui.m.groups))')
    gui.groupListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','Parent', groupBoxV, ...
        'String', gui.Window.UserData.gids,'Value', 1,'Callback', @onGroupListSelection);
    %% MID
    midBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3, 'Visible', 'off');
    uicontrol('Parent', midBoxV,'Style','text','HorizontalAlignment', ...
        'left', 'String','Muscimol Bin:');
    gui.midListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w', ...
        'Parent', midBoxV,'String', {'before','midall','after'}, ...
        'Value', 1,'Callback', @onMidListSelection, 'Enable', 'on');
    if ~gui.m.g(1).after.exists
        set(gui.midListBox,'String', {'before';'midall'});
    end
    %% RUN
    runBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    gui.RunButton = uicontrol( 'Style', 'PushButton','Parent', runBoxV,'String', 'Analysis!', ...
        'Callback', @onRunButton );
    gui.status = uicontrol('Parent', runBoxV,'Style','text', 'Tag', 'status',...
        'HorizontalAlignment', 'left', 'String', 'ok', 'FontSize',9);
    %% PARAMETERS PANEL TIGHTEN UP     
    %Tighten UP
    t = [];
    set(gridBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
    set(gridMidBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
    set(lagBoxV, 'Heights', [20 -1]); t = [t 50];% Make the list fill the space
    set(spikeBinBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
    set(sigmaBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
    set(groupBoxV, 'Heights', [20 -1]); t = [t 200]; % Make the list fill the space
    set(midBoxV, 'Heights', [20 -1]); t = [t 100]; % Make the list fill the space
    set(runBoxV, 'Heights', [30 -1]); t = [t -1]; % Make the list fill the space
    
    set(paramsLayout, 'Heights', t); % Make the lists fill the space
%% CALLBACKS
% GRID
function onGridSliding(hObject,~)
   t = round(get(hObject,'Value')/0.1)*0.1; %nearest 5
   set(gui.gridEdit, 'String', t);
end
function onGridSlide( src, ~ )
    t = round2r(get(src, 'Value'),0.1);
    gui.Window.UserData.gridThresh = t;
    set(gui.gridEdit, 'String', t);
    goodGroups();
    for i = 1:length(gui.Window.UserData.gidsFns)
        f = gui.Window.UserData.gidsFns{i};
        f();
        disp('hrere');
    end
    set(gui.groupListBox,'String', gui.Window.UserData.gids);
    drawnow();
    %run();
end
function onGridEdit( src, ~ )
    t = round2r(str2double(get(src, 'String')),0.1);
    gui.Window.UserData.gridThresh = t;
    set(gui.gridSlider, 'Value', t);
    set(gui.gridEdit, 'String', t);
    goodGroups();
    drawnow();
end
% GRID MID
function onGridMidSliding(hObject,~)
   t = round(get(hObject,'Value')/0.1)*0.1; %nearest 5
   set(gui.gridMidEdit, 'String', t);
end
function onGridMidSlide( src, ~ )
    t = round2r(get(src, 'Value'),0.1);
    gui.Window.UserData.gridThreshMid = t;
    set(gui.gridMidEdit, 'String', t);
    goodGroups();%%%???
    for i = 1:length(gui.Window.UserData.gidsFns)
        f = gui.Window.UserData.gidsFns{i};
        f();
        disp('hrere');
    end
    set(gui.groupListBox,'String', gui.Window.UserData.gids);
    drawnow();
    %run();
end
function onGridMidEdit( src, ~ )
    t = round2r(str2double(get(src, 'String')),0.1);
    gui.Window.UserData.gridThreshMid = t;
    set(gui.gridMidSlider, 'Value', t);
    set(gui.gridMidEdit, 'String', t);
    goodGroups();
    drawnow();
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
    %run();
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
% GROUP
function onGroupListSelection( src, t )
    set(gui.status, 'ForegroundColor',[0,0,1]);
    t = cellstr(get( src, 'String' ));
    gui.m.gid = str2double(t(get( src, 'Value' )));
    gui.m.g = gui.m.groups{gui.m.gid};
    gui.m.g = gui.m.g(gui.Window.UserData.cids{gui.m.gid});
    l = length(gui.m.g);
    set(gui.status, 'String', sprintf('group size: %d',l));
    set(gui.midListBox,'String',{'before';'midall';'after'},'Value', 1);
    if ~gui.m.g(1).after.exists
        set(gui.midListBox,'String', {'before';'midall'});
    end
    %drawnow(); %TRY??
    run();
    %updateInterface();
    %redrawDemo();
end
% MID
function onMidListSelection( src, e )
    %gui.m.mid = get( src, 'Value' )-1;
    t = get( src, 'String');
    gui.m.sesh = str2double(t{get( src, 'Value' )});
    if isnan(gui.m.sesh )
        gui.m.sesh = t{get( src, 'Value' )};
    end
    %run();
    %updateInterface();
    %redrawDemo();
end
% FRAME SLIDER
function onFrameSlide(src, ~)
    gui.m.wini = round(get(src, 'Value'));
    delete(findobj(gui.m.parent,'type','axes'));
    gui.m = driftPlot(gui.m);
    [gui.m.wini, gui.m.title]
    set(gui.status,   'String',  gui.m.title);
    set(gui.viewPanel,'Title',   gui.m.title);
    drawnow();
end
% run
function onRunButton( ~, ~ )
    run();
end % onDemoHelp

%% RUN EXIT
function onExit( ~, ~ )
    % User wants to quit out of the application
    delete( gui.Window );
end % onExit
function run()
    updateM();
    delete(findobj(tab,'type','axes'));
    %m = gui.m; m.fig = gui.Window; g = m.g;  
    set(gui.status, 'String', 'computing....');
    set(gui.status, 'ForegroundColor',[0,0,1]);
    drawnow();
    gui.m.parent = gui.viewContainer;
    plotGroup(gui.m.g, gui.m.parent);
    set(gui.status, 'String', '');
    %enable(handles, true);
end

 %% UPDATE M
function updateM
    gui.m.grid_thresh = gui.Window.UserData.gridThresh;
    gui.m.grid_mid_thresh = gui.Window.UserData.gridThreshMid;
end
%% IDS CIDS
function goodGroups
    
    [gui.Window.UserData.gids, gui.Window.UserData.cids] = groupsByGridThresh(...
        gui.m.groups, gui.Window.UserData.gridThresh, gui.Window.UserData.gridThreshMid);
%     ids = {}; cids = {};
%     for i = 1:length(gui.m.groups)
%         g = gui.m.groups{i};
%         t = [];
%         for j = 1:length(g)
%             if g(j).before.gridscore > gui.Window.UserData.gridThresh...
%             && g(j).midall.gridscore < gui.Window.UserData.gridThreshMid
%                 t(end+1) = j;
%             end
%         end
%         cids{i} = t; 
%         if length(t) >= 2 %only groups with at least 2
%             ids{end+1} = num2str(i);
%         end
%     end
%     gui.Window.UserData.gids = gids;
%     gui.Window.UserData.cids = cids;
end

%run();

%%
end 