function [h,si] = plotRow(data,pd,h,r,si)
    d = plotRowShuffRayData(data,pd);
    [h,si] = plotRowShuffRay(d,pd,h,r,2,si);
    
    
end

function [h,si] = plotRowShuffRay(pdata,pd,h,r,c,si)
    i = 0; ax = []; h = figure(h); %afs = 16; mfs = 16; tfs = 20; lfs = 3;
    afs = 7; mfs = 13; tfs = 9; lfs = 1;
    %both clusters RAYLEIGH
    i=i+1; d = pdata{i}; si=si+1;
    ax(end+1) = subplot(r,c,si); hold on; alpha(ax(end),0.5); %set(ax(end),'color',[.5 .5 .8]);
    plot(ax(end),d.x, d.y,'k.','markersize',mfs); %plot(ax(end),d.x2, d.y2, 'ro','linewidth',1.2);
    %xlim(ax(end),d.xl); ylim(ax(end),d.yl); 
    title(ax(end),d.title,'FontSize',tfs); 
    xlabel(gca,d.xt,'FontSize',afs);ylabel(gca,d.yt,'FontSize',afs)
    %set(ax(end).XTick,'FontSize',afs -2); set(ax(end).yTick,'FontSize',afs -2);%ax(end).YTick.FontSize = afs -2;   
    if isequal(pd.sesh,'PRE')
         %plot(ax(end),d.xdurg, d.ydurg, 'b.','markersize',14); 
         %plot(ax(end),d.xdurl, d.ydurl, 'r.','markersize',14);
    else
         plot(ax(end),[0.4 0.4 1],[1 0.4 0.4],'r--','linewidth',lfs); %plot(ax(end),[0.51 0.51 1],[1 0.51 0.51],'r--');
          text(gca,0.1,1,'Cluster 1','Color','b','FontSize',afs); 
          text(gca,0.6,1,'Cluster 2','Color','r','FontSize',afs);
%         plot(ax(end),d.xdurg, d.ydurg, 'r.','markersize',14); 
%         plot(ax(end),d.xdurl, d.ydurl, 'b.','markersize',14);
    end
    %set(ax(end),'XTickLabel',[0:0.2:1],'fontsize',afs)
    %set(ax(end),'YTickLabel',[0:0.2:1],'fontsize',afs)
    %%%%%%%%%%%%%%%%%%
    
    %less than cluster TEMPORAL
    i=i+1; d = pdata{i}; si=si+1;
    ax(end+1) = subplot(r,c,si); hold on; alpha(ax(end),0.5);
    if isequal(pd.sesh,'PRE')
        p5=plot(ax(end),pd.shufbx,pd.shufmy,'mo','markersize',mfs,'linewidth',lfs);s5 = 'p < %1';
    end
    p1=plot(ax(end),d.x, d.y, 'b.', 'markersize',mfs); 
    title(ax(end),d.title,'FontSize',tfs); xlabel(gca,d.xt,'FontSize',afs);ylabel(gca,d.yt,'FontSize',afs)
    
%     p3=plot(ax(end),[min(d.x) max(d.x)], polyval( polyfit( d.x, d.y,1), [min(d.x) max(d.x)] ),'linewidth',lfs,'color','b');
    
    %text(gca,0.8*max(d.x),polyval(polyfit(d.x, d.y, 1),0.8*max(d.x)),sprintf('%.2f %.2f', round(polyfit(d.x, d.y, 1),2)),'Color','b','FontSize',afs);
    s3 = sprintf('%.1f %.1f', round(polyfit(d.x, d.y, 1),1));
    %greater than cluster
    i=i+1; d2 = pdata{i};  %si=si+1;
    p2=plot(ax(end),d2.x, d2.y, 'r.','markersize',mfs); title(ax(end),d2.title);
%     p4=plot(ax(end),[min(d2.x) max(d2.x)], polyval( polyfit( d2.x, d2.y, 1), [min(d2.x) max(d2.x)] ) ,'linewidth',lfs,'color','r');
    %text(gca,max(d2.x),polyval(polyfit(d2.x, d2.y, 1),max(d2.x)),sprintf('%.2f %.2f', round(polyfit(d2.x, d2.y, 1),2)),'Color','r','FontSize',afs); 
    s4 = sprintf('%.1f %.1f', round(polyfit(d2.x, d2.y, 1),1));
    %xlim(ax(end),[-0.02 0.05]); ylim(ax(end),[-0.02 0.05]);
%     legend(gca,[p1  p3],{'Cluster1',s3},'FontSize',afs);
    %if isequal(pd.sesh,'PRE')
%         legend(gca,[p1 p2 p3 p4 p5],{'Cluster1','Cluster 2',s3,s4,s5},'FontSize',afs);
    %else
        %legend(gca,[p1 p2 p3 p4],{'Cluster1','Cluster2',s3,s4},'FontSize',afs);
    %end
    %xlim(ax(end),d.yl); ylim(ax(end),d.yl);
    %yl = get(ax(end),'yTickLabel');
    %xl = get(ax(end),'XTickLabel');
%     set(ax(end),'XTickLabel',get(ax(end),'XTickLabel'),'fontsize',14)
%     set(ax(end),'YTickLabel',get(ax(end),'XTickLabel'),'fontsize',14)
%     set(ax(end),'XTickLabel',xl,'fontsize',14);
%     set(ax(end),'YTickLabel',yl,'fontsize',14);
    %xlim(ax(end),d.yl); ylim(ax(end),d.yl);
    %set(ax(end),'XTickLabel',[-0.02:0.01:0.03],'fontsize',afs)
    %set(ax(end),'YTickLabel',[-0.02:0.01:0.03],'fontsize',afs)
    for i = 1:length(ax)
           %a = get(ax(i),'XTickLabel');
           
           axis(ax(i),'square'); axis(ax(i),'equal');
    end
    
    

    
    %{
    %hist less than cluster before hist
    i=i+1; d = pdata{i}; si=si+1; 
    ax(end+1) = subplot(r,c,si); hold on; alpha(ax(end),0.2);
    hold on; h1=histogram(ax(end),d.x,'BinEdges',-0.05:0.005:0.05,'FaceColor','b'); %0.15
    xlim(ax(end),d.xl); title(ax(end),d.title); xlabel(gca,d.xt);
    %hist greater than cluster before hist
    i=i+1; d = pdata{i}; %si=si+1;
    %ax(end+1) = subplot(r,c,si); hold on; 
    hold on; h2=histogram(ax(end),d.x,'BinEdges',-0.05:0.005:0.05,'FaceColor','r');  %0.15
    %xlim(ax(end),d.xl); title(ax(end),d.title);
    legend(gca,{'Cluster 1','Cluster 2'});
    %hist less than cluster mid hist
    i=i+1; d = pdata{i}; si=si+1;
    ax(end+1) = subplot(r,c,si); hold on; alpha(ax(end),0.2);
    hold on; histogram(ax(end),d.x,'BinEdges',-0.05:0.005:0.05,'FaceColor','b'); %0.15
    xlim(ax(end),d.xl); title(ax(end),d.title);xlabel(gca,d.xt);
    %hist greater than cluster mid hist
    i=i+1; d = pdata{i}; %si=si+1;
    %ax(end+1) = subplot(r,c,si); hold on; 
    hold on; histogram(ax(end),d.x,'BinEdges',-0.05:0.005:0.05,'FaceColor','r'); %0.15
    xlim(ax(end),d.xl); title(ax(end),d.title);
    legend(gca,{'Cluster 1','Cluster 2'});
    %}
end

function pdata = plotRowShuffRayData(data,pd)
    i1i = 1; i2i = 2; cbi = 3; cmi = 4; iui = 5; r1bi = 6; r2bi = 7; r1mi = 8; r2mi = 9;
    t = data;
    
    if isequal(pd.sesh,'DUR')
        r2 = r2mi;
        r1 = r1mi;
    else
        r2 = r2bi;
        r1 = r1bi;
    end
    
    tg = t(t(:,r2)>=0.4 & t(:,r1)>=0.4,:); tl = t(t(:,r2)<0.4 | t(:,r1)<0.4,:);
    lims = [min(min(t(:,cbi), min(t(:,cmi)))) max(max(t(:,cbi)), max(t(:,cmi)))] ;
    lims = [-0.02 0.03];
    
    ri = 0; pdr = {}; %MID BEF
    
    %1) r1 vs r2 clustering
    r=[];r.stitle=sprintf('%s (%.2f>%.2f) [%s]',pd.sesh, pd.string,pd.gridThreshBef,pd.gridThreshMid);
    r.title=sprintf('%s [%s]','Head Direction of Cell Pairs',pd.sesh);
    r.xt = 'Rayleigh Score Cell 1'; r.yt = 'Rayleigh Score Cell 2';
    r.x = t(:,r1); r.y = t(:,r2); r.co = pd.co; r.xl=[0 1]; r.yl=[0 1]; r.x2 = tg(:,r1); r.y2 = tg(:,r2);
    %%pd.durtg = t(t(:,r2)>=0.5 & t(:,r1)>=0.5,:); pd.durtl = t(t(:,r2)<0.5 | t(:,r1)<0.5,:);
    r.xdurg =  pd.durtg(:,r1); r.ydurg =  pd.durtg(:,r2); r.xdurl =  pd.durtl(:,r1); r.ydurl =  pd.durtl(:,r2);
    
    ri=ri+1; pdr{ri}=r;
    %2) corrs of lower cluster
    r=[];r.title=sprintf('%s','Temporral Correlation Pre vs During');
    r.xt = 'Correlation Cell 1 Cell 2 Pre'; r.yt = 'Correlation Cell 1 Cell 2 During';
    r.x = tl(:,cbi); r.y = tl(:,cmi);
    r.x =  pd.durtl(:,cbi); r.y =  pd.durtl(:,cmi);
    r.co = pd.co; r.xl=lims; r.yl=lims; 
    ri=ri+1; pdr{ri}=r;
    %3) corrs of higher cluster
    r=[];r.title=sprintf('%s','Temporral Correlation Pre vs During');
    r.x = tg(:,cbi); r.y = tg(:,cmi);
    r.x =  pd.durtg(:,cbi); r.y =  pd.durtg(:,cmi); 
    r.co = pd.co; r.xl=lims; r.yl=lims; 
    ri=ri+1; pdr{ri}=r;
    

    
    
    
    
    %4) 5) do hists before mid by cluster (x4) 
    
    r=[];r.title=sprintf('%s','Temporral Correlation Pre Muscimol');
    r.x = tl(:,cbi);r.xl=[-0.05 0.05];  r.xt = 'Corr Cell 1 Cell 2'; %0.15
    ri=ri+1; pdr{ri}=r; 
   
    
    r=[];r.title=sprintf('%s','Temporral Correlation Pre Muscimol');
    r.x = tg(:,cbi);r.xl=[-0.05 0.05]; %0.15
    ri=ri+1; pdr{ri}=r;
    
    r=[];r.title=sprintf('%s','Temporral Correlation During Muscimol');
    r.x = tl(:,cmi);r.xl=[-0.05 0.05]; r.xt = 'Corr Cell 1 Cell 2'; %0.15
    ri=ri+1; pdr{ri}=r;
    
    
    r=[];r.title=sprintf('%s','Temporral Correlation During Muscimol'); 
    r.x = tg(:,cmi);r.xl=[-0.05 0.05]; %0.15
    ri=ri+1; pdr{ri}=r;
    
    pdata = pdr;
end
%{
    t = iud; lims = [min(min(t(:,cbi), min(t(:,cmi)))) max(max(t(:,cbi)), max(t(:,cmi)))] ;
    tg = t(t(:,r2mi)>=0.5 & t(:,r1mi)>=0.5,:); tl = t(t(:,r2mi)<0.5 | t(:,r1mi)<0.5,:);
    
a(end+1) = subplot(row,col,ii);ii=ii+1; plot(t(:,r1mi),t(:,r2mi),'ro'); 
    xlim([0 1]); ylim([0 1]); title(sprintf('decreasing gs pairs (%.2f>%.2f)\nrayleigh i1 vs i2 MID',gridThreshBef,gridThreshMid));  
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cmi),tg(:,cbi),'ro'); title('rayleigh > 0.5 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cmi),tg(:,cbi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    xlim(lims); ylim(lims);
    
a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cmi),tl(:,cbi),'ro'); title('rayleigh < 0.5 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cmi),tl(:,cbi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    xlim(lims); ylim(lims);
    
a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,cbi)); title('hist corrs bef');
    
a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,cmi)); title('hist corrs mid');



    %INTERSECTION UNION
    row = 5; col = 3; ii = 1; 
    figure; 
     t = iud; t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); 
    title(sprintf('decreasing gs pairs (%.2f>%.2f) i/u',gridThreshBef,gridThreshMid)); xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1;plot(tg(:,cmi),tg(:,cbi),'ro');title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cmi),tg(:,cbi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    axis equal;
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cmi),tl(:,cbi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cmi),tl(:,cbi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    axis equal;
    t = iun;t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); title('non decreasing pairs i/u');xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cmi),tg(:,cbi),'ro'); title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cmi),tg(:,cbi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    axis equal;
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cmi),tl(:,cbi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cmi),tl(:,cbi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) );  
    axis equal;
     t = iuo; t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); title('one decreasing pairs i/u');xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cmi),tg(:,cbi),'ro'); title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cmi),tg(:,cbi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    axis equal;
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cmi),tl(:,cbi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cmi),tl(:,cbi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) );
    axis equal;
     t = ium; t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); title('rest of pairs i/u');xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cmi),tg(:,cbi),'ro'); title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cmi),tg(:,cbi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );   
    axis equal;
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cmi),tl(:,cbi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cmi),tl(:,cbi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) );
    axis equal;
     t = iua; t = t(t(:,iui)>0,:); tg = t(t(:,iui)>=0.6,:); tl = t(t(:,iui)<0.6,:); 
    a(end+1) = subplot(row,col,ii);ii=ii+1; hist(t(:,iui),20); title('all pairs i/u');xlim([0 1]);
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cmi),tg(:,cbi),'ro'); title('iu >= 0.6 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cmi),tg(:,cbi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    axis equal;
    a(end+1) = subplot(row,col,ii);         plot(tl(:,cmi),tl(:,cbi),'bo'); title('iu < 0.6 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cmi),tl(:,cbi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) );
    axis equal;


    
    t = findObj(gcf,'type','axes');
    for ri = 1:length(t)
        a = t(ri);
        if mod(ri-1,3) == 0
            mod(ri-1,3)
            axis(a, 'equal');
        end
    end
%}