function TabMisc(fig, tab, m)
    gui.Window = fig;
    gui.m = m;
    
    
    %TOP UI
    AniTabBox = uix.HBoxFlex( 'Parent', tab, 'Spacing', 3);
    controlPanel = uix.Panel('Parent', AniTabBox,'Title', '' );
    gui.viewPanel = uix.Panel('Parent', AniTabBox,'Title', 'Ready','fontweight','bold',...
        'FontSize', 11, 'TitlePosition','centerTop','ForegroundColor','b');
    set( AniTabBox, 'Widths', [200 -1] );
    
    %% VIEW PANEL
    gui.viewContainer = uicontainer('Parent', gui.viewPanel);%'backgroundcolor','k');
    
    %% CONTROL PANEL
    controlBoxV = uix.VBox( 'Parent', controlPanel, 'Spacing', 3);
    t = [];
    %% GROUP
    groupBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 0 ); t = [t -1];
    uicontrol('Style','text','Parent', groupBoxV,'HorizontalAlignment', 'left','String', 'Group:');
    gui.groupListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w','Parent', groupBoxV, ...
        'String', num2str((1:length(gui.m.groups))'),'Value', 1,'Callback', @onGroupListSelection);
    set( groupBoxV, 'Heights', [20 -1] );
    %% MID
    midBoxV = uix.VBox( 'Parent', controlBoxV,'Padding', 3, 'Spacing', 3 ); t = [t 150];
    uicontrol('Parent', midBoxV,'Style','text','HorizontalAlignment', ...
        'left', 'String','Muscimol Bin:');
    gui.midListBox = uicontrol( 'Style', 'list','BackgroundColor', 'w', ...
        'Parent', midBoxV,'String', {'before';'midall';'after'}, ...
        'Value', 1,'Callback', @onMidListSelection, 'Enable', 'on');
    set( midBoxV, 'Heights', [20 -1] );
    set(controlBoxV, 'Heights', t );
    %% CALLBACKS
    %GROUP
    function onGroupListSelection( src, ~ )
        set(gui.viewPanel,'ForegroundColor',[0,0,1], 'Title', '');
        t = cellstr(get( src, 'String' ));
        gui.m.gid = str2double(t(get( src, 'Value' )));
        gui.m.g = gui.m.groups{gui.m.gid};
        l = length(gui.m.g);
        disp(['g: ', num2str(gui.m.gid),' ', num2str(l)]);
        set(gui.viewPanel,'Title', sprintf('group size: %d',l));
        set(gui.midListBox,'String',{'before';'midall';'after'},'Value', 1);
        delete(findobj(tab,'type','axes'));
        set(gui.viewPanel,'Title', sprintf('g: %d length %d', gui.m.gid,l));
        gui.m.parent = gui.viewContainer;
        
        
        %grid scatter
        a = [length(gui.m.g)]; b = [length(gui.m.g)]; c = [length(gui.m.g)];
        for i = 1:length(gui.m.g)
            a(i) = gui.m.g(i).before.gridscore;
            b(i) = gui.m.g(i).midall.gridscore;
            c(i) = gui.m.g(i).after.gridscore;
            %after
        end
        ax = axes('Parent',gui.m.parent);
        plot(ax, 1:length(gui.m.g),[a; b; c],'-.','MarkerSize',32); 
        legend(ax,{'before','muscimol','after'});
        xlabel(ax,' cell'); ylabel(ax,'gridscore');
        title(ax, sprintf('%s','group gridscore [pre, during, post] muscimol '));
        %axis(ax,'equal'); 
        ylim(ax,[-2 2]);  
        set(ax,'ydir','normal','xtick',0:(length(gui.m.g)+1),'xticklabel',0:(length(gui.m.g)+1));
    end
    % MID
    function onMidListSelection( src, ~ )
        t = get( src, 'String');
        disp(t);
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















