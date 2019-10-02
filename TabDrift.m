function TabDrift(fig, tab, m)
    
    gui.Window = fig;
    gui.m = m;
    gui.m.tab = tab;
    gui.m.winsecs = 60;
    gui.m.wini = 1;
    gui.m.sesh = 'before';
    gui.Window.UserData.gidsFns{end+1} = @updateGids;
    
    % Create the panels
    %TOP UI
    DriftTabBox = uix.HBoxFlex( 'Parent', tab, 'Spacing', 3);
    gui.paramPanel = uix.Panel('Parent', DriftTabBox,'Title', 'Parameters' );
    gui.viewPanel = uix.Panel('Parent', DriftTabBox,'Title', 'View','fontweight','bold',...
        'FontSize', 11, 'TitlePosition','centerTop','ForegroundColor','b');
    set( DriftTabBox, 'Widths', [200,-1] );% Adjust the main layout
    %% VIEW PANEL
    viewLayout = uix.VBox( 'Parent', gui.viewPanel,'Padding', 0, 'Spacing', 0);
    t = [];
    gui.viewContainer = uicontainer('Parent', viewLayout); t = [t -1]; %'backgroundColor','w'
    winPanel = uix.Panel('Parent', viewLayout,'Title', 'Window'); t = [t 50];
    %GET MAX MIN
    step = 1;step = step/(10-1); %max-min
    gui.frameSlider = uicontrol( 'Style', 'slider', 'Min', 1,'Max', 10,'Parent', winPanel, ...
        'Value', 1,'SliderStep',[step 2*step],'Callback', @onFrameSlide);
    set(viewLayout, 'Heights', t);
    %% PARAMETERS PANEL   
    paramsLayout = uix.VBox( 'Parent', gui.paramPanel,'Padding', 3, 'Spacing', 3 );
    %WIN
    winBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    uicontrol('Parent', winBoxV,'Style','text',...
        'HorizontalAlignment', 'left','String', 'Window Length(s):');
    winBoxBoxH = uix.HBox( 'Parent', winBoxV,'Padding', 3, 'Spacing', 3 );
    step = 5; wn = 10; wx = 60; step = step/(wx-wn); %max-min
    gui.winSlider = uicontrol( 'Style', 'slider', 'Min', wn,'Max', wx,'Parent', winBoxBoxH, ...
        'Value', wx,'SliderStep',[step 2*step],'Callback', @onWinSlide);
    gui.winEdit = uicontrol( 'Parent', winBoxBoxH, 'Style', 'edit', 'String', wx,...
        'Callback', @onWinEdit);
    addlistener(gui.winSlider,'ContinuousValueChange',@(hObject, event) onWinSliding(hObject, event));
    set( winBoxBoxH, 'Widths', [-4 -1] );   
    %% obsolete LAG
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
    %SPIKE
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
    %SIGMA
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
        'HorizontalAlignment', 'left','String', 'Group:');
    gui.groupListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','Parent', groupBoxV, ...
        'String',gui.Window.UserData.gids,'Value', 1,'Callback', @onGroupListSelection);
    %MID
    midBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    uicontrol('Parent', midBoxV,'Style','text','HorizontalAlignment', ...
        'left', 'String','Muscimol Bin:');
    gui.midListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w', ...
        'Parent', midBoxV,'String', {'before';'midall';'after'}, ...
        'Value', 1,'Callback', @onMidListSelection, 'Enable', 'on');
    if ~gui.m.g(1).after.exists
        set(gui.midListBox,'String', {'before';'midall'});
    end
    %RUN
    runBoxV = uix.VBox( 'Parent', paramsLayout,'Padding', 3, 'Spacing', 3 );
    gui.RunButton = uicontrol( 'Style', 'PushButton','Parent', runBoxV,'String', 'Analysis!', ...
        'Callback', @onRunButton );
    gui.status = uicontrol('Parent', runBoxV,'Style','text', 'Tag', 'status',...
        'HorizontalAlignment', 'left', 'String', 'ok', 'FontSize',9);
    %% PARAMETERS PANEL TIGHTEN UP     
    %Tighten UP
    t = [];
    set(winBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
    set(lagBoxV, 'Heights', [20 -1]); t = [t 50];% Make the list fill the space
    set(spikeBinBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
    set(sigmaBoxV, 'Heights', [20 -1]); t = [t 50]; % Make the list fill the space
    set(groupBoxV, 'Heights', [20 -1]); t = [t 200]; % Make the list fill the space
    set(midBoxV, 'Heights', [20 -1]); t = [t 120]; % Make the list fill the space
    set(runBoxV, 'Heights', [30 -1]); t = [t -1]; % Make the list fill the space
    
    set(paramsLayout, 'Heights', t); % Make the lists fill the space
%% CALLBACKS

function updateGids()
    set(gui.groupListBox, 'String', gui.Window.UserData.gids);
end

% WIN FNS
function onWinSliding(hObject,~)
   t = round(get(hObject,'Value')/5)*5; %nearest 5
   set(gui.winEdit, 'String', t);
end
function onWinSlide( src, ~ )
    t = round2r(get(src, 'Value'),5);
    gui.m.winsecs = t;
    set(gui.winEdit, 'String', t);
    %run();
end
function onWinEdit( src, ~ )
    t = round2r(str2double(get(src, 'String')),5);
    gui.m.windecs = t;
    set(gui.winSlider, 'Value', t);
    set(gui.winEdit, 'String', t);
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
% GROUP FNS
function onGroupListSelection( src, t )
    set(gui.status, 'ForegroundColor',[0,0,1]);
            t = cellstr(get( src, 'String' ));
        gui.m.gid = str2double(t(get( src, 'Value' )));
        gui.m.g = gui.m.groups{gui.m.gid};
    l = length(gui.m.g);
    set(gui.status, 'String', sprintf('group size: %d',l));
    set(gui.midListBox,'String',{'before';'midall';'after'},'Value', 1);
    if ~gui.m.g(1).after.exists
        set(gui.midListBox,'String', {'before';'midall'});
    end    
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

% FRAME SLIDER
function onFrameSlide(src, ~)
    gui.m.wini = round(get(src, 'Value'));
    delete(findobj(gui.m.parent,'type','axes'));
    gui.m = driftPlot(gui.m);
    set(gui.status,   'String',  gui.m.title);
    set(gui.viewPanel,'Title',   gui.m.title);
    drawnow();
end
%
function onRunButton( ~, ~ )
    run();
end % onDemoHelp

%% RUN EXIT
function onExit( ~, ~ )
    % User wants to quit out of the application
    delete( gui.Window );
end % onExit
function run()
    set(gui.status, 'ForegroundColor',[0,0,1]);
    disp(['quest ', num2str(gui.m.gid)]);
    %delete(findobj(gui.Window,'type','axes'));
    delete(findobj(gui.m.tab,'type','axes'));
    drawnow();
    %m = gui.m; m.fig = gui.Window; g = m.g;
    t = get(gui.groupListBox, 'String' );
    gui.m.gid = str2double(t(get( gui.groupListBox, 'Value' )));
    disp(['g: ', num2str(gui.m.gid)])
    gui.m.g = gui.m.groups{gui.m.gid};
    set(gui.status, 'String', 'computing....');
    set(gui.viewPanel,'Title',   'computing....');
    drawnow;
    gui.m.parent = gui.viewContainer;
    gui.m.wini = 1;
    gui.m.windows = {};
    gui.m.grid_thresh = gui.Window.UserData.gridThresh;
    %DRIFT PLOT
    gui.m = driftPlot(gui.m);
    mx = length(gui.m.windows{1}) + 1;
    step = 1; step = step/(mx-1); %max-min
    set(gui.frameSlider,'Min', 1, 'Max', mx, 'Value', 1,'SliderStep',[step 2*step]);
    set(gui.status,   'String',  gui.m.title);
    set(gui.viewPanel,'Title',   gui.m.title);
    drawnow;
    
    %set(gui.status, 'ForegroundColor',[1,0,0]);
    %set(gui.status, 'String', sprintf('good cells[%s,%d]','',length(gui.m.g)));
    %set(gui.viewPanel,'Title',sprintf('good cells[%s,%d]','',length(gui.m.g)));
    %{
    if isempty(gui.m.good) || length(gui.m.good) < 2
        set(gui.status, 'ForegroundColor',[1,0,0]);
        set(gui.status, 'String', sprintf('good cells[%s,%d]',v,length(good)));
        set(gui.viewPanel,'Title',sprintf('good cells[%s,%d]',v,length(good)));
    else
        mx = length(gui.m.windows{1}) + 1;
        step = 1; step = step/(mx-1); %max-min
        set(gui.frameSlider,'Min', 1, 'Max', mx, 'Value', 1,'SliderStep',[step 2*step]);
        set(gui.status,   'String',  gui.m.title);
        set(gui.viewPanel,'Title',   gui.m.title);
    end
    %}
end

%run();

%%
end 