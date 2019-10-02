% uncomment to load cells
%load(sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45));
%load('Z:\\data\\noam\\Data\\muscimol\\noam\\cells_45min_d_patchtraj_rayleigh'); %lab
%load('C:\Noam\Data\muscimol\cells15nan');
%add more to quantify gridscpre ??
%preamble
function f2(cellsn)
dbstop if error    

ci=[171,174,176,177]; 
c1 = cellsn(ci(1)); c2 = cellsn(ci(3)); %171 176 %  c1 34; c2 35; 112 119
%ci = [112,114,117,119]; %[171,174,176,177];    %ci = [34, 35, 37, 41];,113 116,118

fig = figure(992); clf; tic
set(fig,'Color','w', 'Position', [200 50 800 800]); p = {}; s1 = {};
gtop = uix.GridFlex('Parent', fig,'Spacing',5, 'BackgroundColor','w','DividerMarkings','on');

%% A
gA = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axA = {}; fs = 14; p = [];
axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
%Text
text(0,0.5,'A','fontweight','bold','fontsize',fs);
axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
text(0.5,0.5,'Pre', 'fontsize',fs-2);
axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'),'visible','off');
text(0.5,0.5,'During','fontsize',fs-2);
%TRAJ 1
uix.Empty('Parent',gA);

s1.x=''; s1.y = s1.x; s1.t = '';
axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w')); 
axA{end} = plotTR(axA{end},c1.before,s1);    
axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w'));
axA{end} =plotTR(axA{end},c1.midall,s1);
%TRAJ 2
uix.Empty('Parent',gA);
axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w')); 
axA{end} = plotTR(axA{end},c2.before,s1);    
axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w')); 
axA{end} =plotTR(axA{end},c2.midall,s1);
%TCC
p.lag = 2000; p.movmean = 100;  
uix.Empty('Parent',gA);
axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w')); 
axA{end}  = plotTimeCorr(axA{end},c1.before,c2.before,p,s1);
axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w')); 
axA{end}  = plotTimeCorr(axA{end},c1.midall,c2.midall,p,s1);
%SCC
uix.Empty('Parent',gA);
axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w')); 
axA{end}  = plotSpaceCorr(axA{end},c1.before,c2.before,p,s1);
axA{end+1} = axes('Parent',uicontainer('Parent',gA,'BackgroundColor','w')); 
axA{end}  = plotSpaceCorr(axA{end},c1.midall,c2.midall,p,s1);

uix.Empty('Parent',gA);

set(gA,'Heights', [25 -1 -1]);%,'Widths', [-1 -1 -1]);
 


%% C SAME GROUP
gB = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
axB = {};  p.movmean = 100;  
axB{end+1} = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'),'visible','off');
text(0,0.5,'B','fontweight','bold','fontsize',fs);
s1.x=''; s1.y = ''; s1.t = ''; %p.off = 0;
for i = 1:len(ci)-1
    for ii = i+1:len(ci)
        c1 = cellsn(ci(i));c2 = cellsn(ci(ii));
        if i+ii>3
            uix.Empty('Parent',gB);
        end
        axB{end+1} = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'));
        axB{end}  = plotTimeCorr(axB{end},c1.before,c2.before,p,s1);
        axB{end+1} = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w'));
        axB{end}  = plotTimeCorr(axB{end},c1.midall,c2.midall,p,s1); 

        c1 = cellsn(ci(i));c2 = cellsn(ci(ii));
        axB{end+1} = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w')); 
        axB{end}  = plotSpaceCorr(axB{end},c1.before,c2.before,p,s1);
        axB{end+1} = axes('Parent',uicontainer('Parent',gB,'BackgroundColor','w')); 
        axB{end}  = plotSpaceCorr(axB{end},c1.midall,c2.midall,p,s1);
        
    end
end
set(gB,'Heights', [25 -1 -1 -1 -1]);%,'Widths', [-1 -1 -1]);


%% epilogue
% set(gtop,'Heights', [-1.6 -0.97 -3]);
set(gtop,'Heights', [-1 -2]);

end



%% ?
%% B SMOOTHING
%     gC = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
%     c1 = cellsn(34); c2 = cellsn(35); axC = {}; 
%     axC{end+1} = axes('Parent',uicontainer('Parent',gC,'BackgroundColor','w'),'visible','off');
%     text(0,0.5,'B','fontweight','bold','fontsize',fs);
% %    	p.lag = 5*1000;
%     mm = [1, 10, 50, 250, 1000, 2000];
%     bc = []; mc= [];
%     for i = 1:len(mm)
%         if i>1
%             uix.Empty('Parent',gC);
%         end
%         p.movmean = mm(i);s1.t = sprintf('%dms',mm(i));
%         axC{end+1} = axes('Parent',uicontainer('Parent',gC,'BackgroundColor','w')); 
%         axC{end}  = plotTimeCorr(axC{end},c1.before,c2.before,p,s1);
%     end
%     set(gC,'Heights', [25 -1]);%,'Widths', [-1 -1 -1]);

