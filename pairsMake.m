function [pairs,cels,group,chd,phd,chdn,phdn]=pairsMake(cellsn, bthresh, mthresh)

%load
%load('Z:\data\noam\muscimol\cells15nan');
%load('C:\Noam\Data\muscimol\cells15nan');
%load('.\Data\gshuff 06-Aug-2019-547 with aggr');
sss='pearson with nans';
%bthresh=0.5; mthresh=0.2; %V2
b=[cellsn.before]; l =     arrayfun(@(z) len(z.st)>=100,b)'; 
m=[cellsn.midall]; l = l & arrayfun(@(z) len(z.st)>=100,m)'; 
gss='gridscore'; acss='ac';
b=[b.(gss)]'; m=[m.(gss)]';
gth= b>bthresh & m<mthresh; %gth= b>0 & m<0; 
%ppgbma = pgbma <= ceil(max(pgbma(:))*0.01); %%<<%% picking percentile significance
gsig=ones(len(cellsn),1);%gsig= (ppgbma(:,1));% & ~ppgbma(:,2)); %significant gridscores
gclls=cellsn(gsig&l&gth); clear gth; clear l;clear b; clear m;
%sss=[sss ' | selection: pre pos sign gscore | dur neg gscore |'];
sss=[sss ' | selection: pre ' n2(bthresh) ' | dur ' n2(mthresh) ' |'];

%group by recording date and animal
gs =                 unique(str2double(strcat({gclls.id},{gclls.date}))    )';
group = arrayfun(@(x) gclls(str2double(strcat({gclls.id},{gclls.date}))==x ),gs,'uni',0);
groupis= cellfun (@(x) [x.ind],group,'uni',0)';
%sort by group id;
[~,i]=sort(cellfun(@(x) x(1),groupis));
group=groupis(i)';
clear groupis; clear gs; clear i; clear gccls;

threshRemoveSpikes=-1; 
threshRemoveCell=0.05;
percentremoved=[];

for gi = 1:len(group)
    gcis=group{gi};
    g = cellsn(gcis);
    rcs = [];
    for j = 1:len(g)-1
        for k = j+1:len(g)
            c1 = g(j).before; c2 = g(k).before;
            s1=floor(c1.st*1000)+1; s2=floor(c2.st*1000)+1; min_time=min(min(s1),min(s2)); %TIME IN MS
            offset = 0; s1 = s1-min_time+offset; s2 = s2-min_time+offset; max_time=max(max(s1),max(s2))+offset;
            [s1i, s2i] = removeOverlappingSpikes(s1,s2, 1); rs1=s1(s1i); rs2=s2(s2i);
            %fprintf('####b.g%2d i%3dxi%3d: %3d :: %3d\n',ri, g(j).ind ,g(k).ind,length(s1)-length(rs1),length(s2)-length(rs2));
            rem = length(s1)-length(rs1); %assert(abs( (length(s1)-length(rs1)) - (length(s2)-length(rs2)) )<3)
            brem = max(rem/length(s1),rem/length(s2));
            if brem>threshRemoveCell || false
                fprintf('b.g%2d: spikes removed %3dx%3d = %3d, --- c1 %2d%% c2 %2d%% \n',gi, g(j).ind,g(k).ind,rem,...
                    100*round(1-length(rs1)/length(s1),2),100*round(1-length(rs2)/length(s2),2));
            else
                percentremoved(end+1)=brem;
            end
            if brem >  threshRemoveCell
                %remove cell with higher percentage spikes removed
                rid = gcis(j);
                %s2 got reduced more
                if length(rs2)/length(s2) < length(rs1)/length(s1)
                    rid = gcis(k);
                end
                %fprintf('g.%2d: b c.%3d: %3f%%\n',gids(ri), rid, 100*max(brem,mrem));
                rcs = [rcs rid];
            end
            c1 = g(j).midall; c2 = g(k).midall;
            s1=floor(c1.st*1000)+1; s2=floor(c2.st*1000)+1; min_time=min(min(s1),min(s2)); %TIME IN MS
            offset = 0; s1 = s1-min_time+offset; s2 = s2-min_time+offset; max_time=max(max(s1),max(s2))+offset;
            [s1i, s2i] = removeOverlappingSpikes(s1,s2, 1); rs1=s1(s1i); rs2=s2(s2i);
            rem = length(s1)-length(rs1); assert(abs( (length(s1)-length(rs1)) - (length(s2)-length(rs2)) )<=3)
            mrem = max( rem/length(s1), rem/length(s2) );
            if mrem>threshRemoveCell || false
                fprintf('m.g%2d: spikes removed %3dx%3d = %3d, --- c1 %2d%% c2 %2d%% \n',gi, g(j).ind,g(k).ind,rem,...
                        100*round(rem/length(s1),2),100*round(rem/length(s2),2));        
                       %100*round(1-length(rs1)/length(s1),2),100*round(1-length(rs2)/length(s2),2));
            else
                percentremoved(end+1)=mrem;
            end
            if mrem > threshRemoveCell
                %remove cell with higher percentage spikes removed
                rid = gcis(j);
                %s2 got reduced more
                if length(rs2)/length(s2) < length(rs1)/length(s1)
                    rid = gcis(k);
                end
                %fprintf('g.%2d: m c.%3d: %3f%%\n',gids(ri), rid, 100*max(brem,mrem));
                rcs = [rcs rid];
            end
        end
    end
    %remove cells
    if ~isempty(rcs)
        rcs = unique(rcs);
        fprintf('$ removing %d cells from g.%d\n',length(rcs),gi);
        [rcs]
        group{gi} = setdiff(gcis, rcs);
    end
end
clear gi;clear gcis;clear g;clear rcs;clear rid;clear rem;clear brem;clear mrem;clear offset;
clear c1;clear c2;clear s1;clear s2;clear s1i;clear s2i;clear rs1;clear rs2;clear max_time;
clear min_time; clear threshRemoveSpikes;clear threshRemoveCell;


%remove 1 cell groups
group=group(cellfun(@(x) len(x)>1,group));
%get pairs
pairs=cellfun(@(x) nchoosek(x,2),group,'uni',false);
pairs=vertcat(pairs{:});
cels=unique(pairs(:));

[chd, phd, chdn,phdn ] =clusters(0.4,0.4,pairs,cellsn,cels); %SAME AS f4

['cells pairs hdcells hdpairs']
[len(cels) len(pairs) len(chd) len(phd)]
sss=[sss ' cells ' n2(len(cels)) ' pairs ' n2(len(pairs)) ' hd cells ' n2(len(chd)) ' hd pairs ' n2(len(phd))];

%%{
%CALLING PLOTS
%fp(1,sss,pairs,cellsn,acss,gss); fp(0,sss,pairs,cellsn,acss,gss);

%}
clear acss;clear gclls;clear gq;clear gr;clear gsig;clear gss;clear gt;clear j;clear nl;clear k;clear p;
clear pgb;clear pgm;clear ax;clear s1;clear fig;%clear ;clear ;


function [chd, phd, chdn,phdn ]=clusters(rsth,rstl,pairs,cellsn,cels)
%   rsth=0.4; rstl=0.4;    
rsp=[];rs=[]; 
for i=1:length(pairs)
    rsp(i,:) = [cellsn(pairs(i,1)).before.rayleigh_score,...
                cellsn(pairs(i,2)).before.rayleigh_score,...
                cellsn(pairs(i,1)).midall.rayleigh_score,...
                cellsn(pairs(i,2)).midall.rayleigh_score];
end
for i=1:length(cels)
    rs(i,:)  = [cellsn(cels(i)).before.rayleigh_score,...
        cellsn(cels(i)).midall.rayleigh_score];
end
ccl2 = rs(:,2)>=rsth & rs(:,1)<rstl; 
ccl1 = rs(:,2) <rstl & rs(:,1)<rstl; %ccl1 = ~ccl2; 
%rsp: rscore[c1b c2b c1m c2m]
cl2 = rsp(:,3)>=rsth & rsp(:,4)>=rsth ...   %midall
    & rsp(:,1)< rstl & rsp(:,2)< rstl;      %before
cl1 = rsp(:,3)< rstl & rsp(:,4)< rstl ...     %midall
    & rsp(:,1)< rstl & rsp(:,2)< rstl;     %before
%cl1 = ~cl2;
chd = cels(ccl2);
chdn = cels(ccl1);
i=1:len(pairs);
phd =  i(cl2)';%pairs(cl2,:);
phdn = i(cl1)'; 


clear rs;clear rsp;clear rsth;clear rstl;clear ccl1;clear ccl2; clear cl1;clear cl2; 
clear i; clear rsth; clear rstl;

end

function fp(sh,sss,pairs,cellsn,acss,gss)
cels=unique(pairs(:));
if sh;ss='bef ';fig = figure(21);set(fig,'Color','w', 'Position', [40     50 900 900]);
else; ss='mid ';fig = figure(22);set(fig,'Color','w', 'Position', [50+900 50 900 900]);end;clf
gt = uix.GridFlex('Parent', fig,'Spacing',5, 'BackgroundColor','w','DividerMarkings','on');
gq = uix.GridFlex('Parent', gt ,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
ax=axes('Parent',uicontainer('Parent',gq,'BackgroundColor','w'),'visible','off');
ss=[ss sss];
text(ax,0,0.5,ss,'fontweight','bold','fontsize',12);
gr = uix.GridFlex('Parent', gt,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
for i = toRow(cels)
    c=cellsn(i);
axes('Parent',uicontainer('Parent',gr,'BackgroundColor','w'),'visible','off');
if sh 
    imgsc(c.before.(acss),2);axis off;title(['i' n2(c.ind) 'pre' n2(c.before.(gss),1)]);
else
    imgsc(c.midall.(acss),2);axis off;title(['i' n2(c.ind) 'mid' n2(c.midall.(gss),1)],'color','r');
end
end
nl = ceil(sqrt(len(cels)));if mod(nl, 2)==1;nl=nl-1;end
set(gr,'Heights',zeros(1,nl)-1);
set(gt,'Heights',[25 -1]);
end

end

% %mid
% fig = figure(22); clf; gr=[]; gt=[];
% set(fig,'Color','w', 'Position', [20 60 1800 900]); p = {}; s1 = {};
% gt = uix.GridFlex('Parent', fig,'Spacing',5, 'BackgroundColor','w','DividerMarkings','on');
% gq = uix.GridFlex('Parent', gt ,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
% ax=axes('Parent',uicontainer('Parent',gq,'BackgroundColor','w'),'visible','off');
% sss=[sss ' cells ' n2(len(cels)) ' pairs ' n2(len(pairs)) ' hd cells ' n2(len(chd)) ' hd pairs ' n2(len(phd))];
% text(ax,0,0.5,['MID ' sss],'fontweight','bold','fontsize',12);
% gr = uix.GridFlex('Parent', gt,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
% for i = toRow(cels)
%     c=cellsn(i);
% %axes('Parent',uicontainer('Parent',gr,'BackgroundColor','w'),'visible','off');
% %imgsc(c.before.(acss),2);axis off;title(['i' n2(c.ind) 'pre' n2(c.before.(gss),1)]); 
% axes('Parent',uicontainer('Parent',gr,'BackgroundColor','w'),'visible','off'); 
% imgsc(c.midall.(acss),2);axis off;title(['i' n2(c.ind) 'mid' n2(c.midall.(gss),1)],'color','r');
% end
% nl = ceil(sqrt(len(cels)));if mod(nl, 2)==1;nl=nl-1;end
% set(gr,'Heights',zeros(1,nl)-1);
% set(gt,'Heights',[25 -1]);

%end
%{
figure(834)
 for i=7
   i
   c=cellsn(i).before;
   subplot(222)
   imgsc(c.rm,1); title(len(c.st));
   subplot(221)
   imgsc(c.ac,2); title(gridscore2(c.ac,2));
   c=cellsn(i).midall;
   subplot(223)
   imgsc(c.ac,2); title(['midall ' n2(c.gridscore)]);
   subplot(224)
   imgsc(c.rm,1); title(['midall shi' n2(pgbma(cellsn(i).ind,2))]);
   
   suptitle(n2(cellsn(i).ind))   
end

%}


% for i=1:len(gclls)
%    i
%    a=gclls(i).before;
%    subplot(222)
%    imgsc(a.rm,1); title(len(a.st));
%    subplot(221)
%    imgsc(a.ac2,2); title(gridscore2(a.ac2,2));
%    a=gclls(i).midall;
%    subplot(223)
%    imgsc(a.ac2,2); title(['midall ' n2(a.gs2)]);
%    subplot(224)
%    imgsc(a.rm,1); title(['midall shi' n2(pgbma(gclls(i).ind,2))]);
%    
%    suptitle(n2(gclls(i).ind))
%    
% end
