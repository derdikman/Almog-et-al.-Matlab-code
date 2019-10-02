function GuiDamo()

%GUI = Figure + data structure

gui = struct();

gui.window = figure(...
    'Position', [50,50, 800, 400],'Name', 'gui time !', ...
    'NumberTitle', 'off','MenuBar', 'none','Toolbar', 'none');

%global data
gui.window.UserData.global = 42;

tabs = uix.TabPanel( 'Parent',gui.window);

gui.plotTab  = uix.Panel('Parent', tabs);
gui.aniTab   = uix.Panel('Parent', tabs);

tabs.TabTitles = {'plot','replay'};

tabs.Selection = 2; % <<<< START UP TAB                

plotTab(gui.window, gui.plotTab);
aniTab(gui.window, gui.aniTab);

end

function plotTab(fig, tab)

tabdata = struct();
tabdata.t = 0:0.01:(2*pi);
tabdata.x = 1;

tabHBox = uix.HBoxFlex( 'Parent', tab, 'Padding', 2, 'Spacing',3);
controlPanel = uix.Panel('Parent', tabHBox,'Title', 'parameters');
%,'FontSize', 11, 'TitlePosition','centerTop');
viewPanel = uix.Panel('Parent', tabHBox,'Title', 'the view panel');
set( tabHBox, 'Widths', [200 -1] );

%% CONTROL PANEL
controlBoxV = uix.VBox( 'Parent', controlPanel, 'Spacing', 3);
%slider
slider = uicontrol(  'Parent',controlBoxV,'Style', 'slider', 'Min', 1,'Max', 10, 'Value',tabdata.x,...
    'Callback', @onSlide); %'SliderStep',[1/9 3/9],
%run button
runButton = uicontrol('Parent', controlBoxV,'Style', 'PushButton','String', 'Yala','Callback',@onRunButton);
%set heights of ui elements
set(controlBoxV,'Heights',[20,50]);


%% CALLBACKS
    function onSlide(src,event)
        tabdata.x = get(src, 'Value')
        run()
    end

    function onRunButton( ~, ~ )
        run();
    end % onDemoHelp

    function run()
        delete(findobj(tab,'type','axes'));
        ax = axes('Parent',viewPanel);
        plot(ax,tabdata.t, sin(tabdata.x*tabdata.t));
        title(ax,sprintf('x = %.1f',tabdata.x));
        axis(ax,'square');axis(ax,'tight');
        
    end

end

function aniTab(fig, tab);
tabdata = struct();
tabdata.t0 = 0:0.01:(2*pi);
tabdata.t = tabdata.t0;
tabdata.u = 1;
tabdata.v = 1;
tabdata.w = 1;
tabdata.x = 1;
tabdata.y = 1;
tabdata.z = 1;
tabdata.a = 0;
tabdata.stop = false;

tabHBox = uix.HBoxFlex( 'Parent', tab, 'Padding', 2, 'Spacing',3);
controlPanel = uix.Panel('Parent', tabHBox,'Title', 'parameters');
%,'FontSize', 11, 'TitlePosition','centerTop');
viewPanel = uix.Panel('Parent', tabHBox,'Title', 'the view panel');
set( tabHBox, 'Widths', [200 -1] );

%% CONTROL PANEL
controlBoxV = uix.VBox( 'Parent', controlPanel, 'Spacing', 3);
%slider
sh = 20; t = [];
uicontrol('Parent',controlBoxV,'Style','slider','Min',1,'Max',100,'Value',tabdata.u,'Callback',@onUSlide);t=[t sh]; %'SliderStep',[1/9 3/9],
uicontrol('Parent',controlBoxV,'Style','slider','Min',1,'Max',100,'Value',tabdata.v,'Callback',@onVSlide);t=[t sh];
uicontrol('Parent',controlBoxV,'Style','slider','Min',1,'Max',100,'Value',tabdata.w,'Callback',@onWSlide);t=[t sh];
uicontrol('Parent',controlBoxV,'Style','slider','Min',1,'Max',100,'Value',tabdata.x,'Callback',@onXSlide);t=[t sh];
uicontrol('Parent',controlBoxV,'Style','slider','Min',1,'Max',100,'Value',tabdata.y,'Callback',@onYSlide);t=[t sh];
uicontrol('Parent',controlBoxV,'Style','slider','Min',1,'Max',100,'Value',tabdata.z,'Callback',@onZSlide);t=[t sh];
uicontrol('Parent',controlBoxV,'Style','slider','Min',-3,'Max',3,'Value', tabdata.a,'Callback',@onASlide);t=[t sh];
%run button
uicontrol('Parent', controlBoxV,'Style', 'PushButton','String', 'Yala','Callback',@onRunButton);
uicontrol('Parent', controlBoxV,'Style', 'PushButton','String', 'Stop!','Callback',@onStopButton);
%set heights of ui elements
set(controlBoxV,'Heights',[t,-1,-1]);
    function onUSlide(src,event)
        tabdata.u = get(src, 'Value')
    end
    function onVSlide(src,event)
        tabdata.v = get(src, 'Value')
    end
    function onWSlide(src,event)
        tabdata.w = get(src, 'Value')
    end
    function onXSlide(src,event)
        tabdata.x = get(src, 'Value')
    end
    function onYSlide(src,event)
        tabdata.y = get(src, 'Value')
    end
    function onZSlide(src,event)
        tabdata.z = get(src, 'Value')
    end
    function onASlide(src,event)
        tabdata.a = get(src, 'Value')
        tabdata.t = tabdata.t0*10^tabdata.a;
    end
    function onRunButton( ~, ~ )
        run();
    end % onDemoHelp
    function onStopButton( ~, ~ )
        tabdata.stop = true;
    end % onDemoHelp
    
    function run(~,~)
        delete(findobj(tab,'type','axes'));
        ax = axes('Parent',viewPanel);
        tabdata.stop = false;
        %title(ax,sprintf('x = %.1f',tabdata.x));
        i = 1;
        while true
            i = i + 1;
            if(i>length(tabdata.t))
                i = 1;
            end
            t = [tabdata.t(i+1:length(tabdata.t)), tabdata.t(1:i)];
            a = sin(tabdata.u.*t).*sin(tabdata.x.*t)*100;
            b = sin((tabdata.y/2).^(tabdata.v.*t))*100;
            c = sin(tabdata.w.*t).*log(t*tabdata.z)*100;
            plot(ax,tabdata.t',...       
                [a+b+c...
                %sin(tabdata.u.*t);...
                %sin(tabdata.v.*t);...
                %sin(tabdata.w.*t);...
                %sin(tabdata.x.*t);...
                %sin(tabdata.y.*t);...
                %sin(tabdata.z.*t)
                ]');
            %axis(ax,'square');
            axis(ax,'tight');
            sound(a+b+c);
            
            %sound(sin(tabdata.x.*t)*100);
            %sound(sin(tabdata.y.*t)*100);
            %sound(sin(tabdata.z.*t)*100);
            drawnow;
            if tabdata.stop
                tabdata.stop = false;
                delete(findobj(tab,'type','axes'));
                ax = axes('Parent',viewPanel);
                break
            end
        end
    end

end
