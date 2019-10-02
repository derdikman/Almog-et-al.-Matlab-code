%{

%}

function   shufflingMain()
    warning off;
    %%%
    binsize = 45;
    %fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_c_midall_gridscore.mat',binsize); first set
    fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45);
    load(fn);
    %PRINT OUT AREAS OF RECORDING
    for i = 1:length(cells)
        t = cells{i};
        fprintf('%d: [%s%s] -  %s\n',t.ind,t.id,t.date,t.area);
    end  
    iCbmIuRbm250u = {};
    [groups ~] = findSimultaneouslyRecordedCells(cells);
    gridThreshBef = 0.3; gridThreshMid = 0.25;        
    [gids, cids] = groupsByGridThresh(groups, gridThreshBef, gridThreshMid);
    %get group mean gridscore drop
    gsd = []; gsdd = []; ac = {}; %SUMMARY STATS
    for gi = 1:length(gids)
        g = groups{gids(gi)}; g = g(cids{gids(gi)});
        for ci = 1:length(g)
            ac{end+1} = struct2cell(g(ci));
            %{
            gsd(end+1) = g(ci).before.gridscore - g(ci).midall.gridscore;
            [g(ci).before.gridscore g(ci).midall.gridscore]
            if g(ci).midall.gridscore ~= -2
                gsdd(end+1) = g(ci).before.gridscore - g(ci).midall.gridscore;
            end
            %}
        end
    end
    %SUMMARY STATS
    av = zeros(length(ac),7);
    for i = 1:length(ac)
       t = ac{i}
       av(i,:) = [t{1} t{2}.max_r t{3}.max_r t{2}.gridscore t{3}.gridscore t{2}.rayleigh_score t{3}.rayleigh_score]; 
    end
    av(av(:,5)==-2,5) = 0;
    subplot(321); histogram(av(:,2),'BinLimits',[0,200],'BinWidth',10); title('Pre');    xlabel('max firing rate(Hz)'); 
    subplot(322); histogram(av(:,3),'BinLimits',[0,200],'BinWidth',10); title('During'); xlabel('max firing rate(Hz)');
    subplot(323); histogram(av(:,4),'BinLimits',[-2,2],'BinWidth',0.2); xlabel('grid score'); 
    subplot(324); histogram(av(:,5),'BinLimits',[-2,2],'BinWidth',0.2); xlabel('grid score'); 
    subplot(325); histogram(av(:,6),'BinLimits',[0,1],'BinWidth',0.05); xlabel('rayleigh score'); 
    subplot(326); histogram(av(:,7),'BinLimits',[0,1],'BinWidth',0.05); xlabel('rayleigh score');
    suptitle('Decreasing Grid Score Cells n=69');
   
    
     for ri = 1:length(gids)
        gis = gids(ri);
        cis = cids{gis};
        g = groups{gis};
        g = g(cis);
        tic
        for j = 1:length(g)-1
            for k = j+1:length(g)
                cb = corrBefore{g(j).ind,g(k).ind};cb=cb(1);
                cm = corrMid{g(j).ind,g(k).ind};cm=cm(1);
                %cb = shuffleTimeCorrelations (g.before,g(k).before,p);
                %cm = shuffleTimeCorrelations (g(j).midall,g(k).midall,p);
%                 iu = Intersection_Over_Union(g(j).before.module.x, g(j).before.module.y, g(k).before.module.x,g(k).before.module.y);
                iCbmIuRbm250u{g(j).ind,g(k).ind} = [g(j).ind g(k).ind cb cm -1 ... 
                    g(j).before.rayleigh_score g(k).before.rayleigh_score g(j).midall.rayleigh_score g(k).midall.rayleigh_score];
            end
        end
    end
     i1i = 1; i2i = 2; cbi = 3; cmi = 4; iui = 5; r1bi = 6; r2bi = 7; r1mi = 8; r2mi = 9;
     
     iud = []; iun = []; iuo = []; ium = []; iua = []; iug = []; iudd = [];                                       
     for ri = 1:length(groups)
       gis = gids(ri);
        cis = cids{gis};
        g = groups{gis};
        g = g(cis);
        for j = 1:length(g)-1
            for k = j+1:length(g)
                i1 = g(j).ind;  i2 = g(k).ind; gm1 =  g(j).midall.gridscore; gm2 =  g(k).midall.gridscore;
                c1 = g(j).before; c2 = g(k).before; cb = iCbmIuRbm250u{i1,i2}(3); cm = iCbmIuRbm250u{i1,i2}(4);
                iu = Intersection_Over_Union(c1.module.x,c1.module.y, c2.module.x,c2.module.y);
                %t = [iu i1 i2 cb cm];
                t = iCbmIuRbm250u{i1,i2};
                %both cells are grid
                if c1.gridscore >= gridThreshBef  && c2.gridscore >= gridThreshBef
                        iug(end +1,:) = t;
                        %both decreasing
                    if     gm1 <= gridThreshMid  && gm2 <= gridThreshMid
                        iud(end +1,:) = t;
                        %both non decreasing
                    elseif gm1 >  gridThreshMid && gm2 >  gridThreshMid
                        iun(end +1,:) = t;
                        %one decreasing
                    elseif(gm1 <= gridThreshMid && gm2 >  gridThreshMid) ||...
                          (gm1 >  gridThreshMid && gm2 <  gridThreshMid)
                        iuo(end +1,:) = t;
                        %should be empty
                    else
                        ium(end +1,:) = t;
                    end
                end
                iua(end +1,:) = t;
            end
        end
     end
     iudd = [iun; iuo]; %grid but not decreasing
    %iCbmIuRbm250u              i1 i2 cb cm iu r1b r2b r1m r2m 
                                 %1  2  3  4  5  6   7   8   9                                   

    %RAYLEIGH CLUSTER PLOT ROW
    pd.gridThreshBef = gridThreshBef; pd.gridThreshMid = gridThreshMid;
    pd.durtg = iud(iud(:,r2mi)>=0.4 & iud(:,r1mi)>=0.4,:); pd.durtl = iud(iud(:,r2mi)<0.4 | iud(:,r1mi)<0.4,:);
    
    h=figure(1);clf(h);set(h,'color','w'); si = 0; row = 2;  %h=figure('position',[-1700,100,1600,800]); 
    pd.sesh = 'DUR'; pd.co = 'ro'; pd.string = 'decreasing gs pairs'; [h,si] = plotRow(iud,pd,h,row,si);
    pd.sesh = 'PRE'; pd.co = 'bo'; pd.string = 'decreasing gs pairs'; [h,si] = plotRow(iud,pd,h,row,si);
    
    
    
    ax = axes('Parent',h,'position',[0.5,0.98,0,0]);axis(ax,'off');
    pt = {'FontSize',12, 'fontweight','bold', 'horizontalalignment','center'}; 
    %text(ax,0,0,'Head Direction and Temporal Correlation of Simultaneously Recorded Grid Cells',pt{:});    
    %pd.sesh = 'MID'; pd.co = 'ro'; pd.string = 'non decreasing gs pairs';[h,si] = plotRow(iudd,pd,h,row,si);
    %pd.sesh = 'BEF'; pd.co = 'bo'; pd.string = 'decreasing gs pairs'; [h,si] = plotRow(iud,pd,h,row,si);
    %pd.sesh = 'BEF'; pd.co = 'bo'; pd.string = 'non decreasing gs pairs';[h,si] = plotRow(iudd,pd,h,row,si);
    
%%%%%%%%%%%%%
%%%%%%%%%%%%%%
    
    naxes = 6; ncol = 2;
    nRow = ceil( naxes / ncol ) ;
    rowH = 0.75 / nRow ;  colW = 0.75 / ncol ; %width of plot
           %offset
    colX = 0.05 + linspace( 0, 0.9, ncol+1 ) ; colX = colX(1:end-1) ;
    rowY = 0.05 + linspace( 0.9, 0, nRow+1 ) ; rowY = rowY(2:end);    
    %rowId = ceil( ai / ncol );
    %colId = ai - (rowId - 1) * ncol;
    %x = colX(colId); y = rowY(rowId);
    %p = [x, y, colW, rowH];
     x = colX( ai - (ceil( ai / ncol ) - 1) * ncol); y = rowY(ceil( ai / ncol ));p = [x, y, colW, rowH];
    
    
    %RAYLEIGH 
    row = 4; col = 2; ii =1; figure; a = [];
    t = iud; lims = [min(min(t(:,cbi), min(t(:,cmi)))) max(max(t(:,cbi)), max(t(:,cmi)))] ;
    tg = t(t(:,r2mi)>=0.5 & t(:,r1mi)>=0.5,:); tl = t(t(:,r2mi)<0.5 | t(:,r1mi)<0.5,:);  
    %r1 vs r2
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(t(:,r1mi),t(:,r2mi),'go'); 
    xlim([0 1]); ylim([0 1]); title(sprintf('decreasing gs pairs (%.2f>%.2f)\nrayleigh i1 vs i2 MID',gridThreshBef,gridThreshMid));  
    ii = ii+1;
    %cb vs cr less than cluster
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tl(:,cbi),tl(:,cmi),'ro'); title('rayleigh < 0.5 corr b vs m');
    hold on; plot([min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    xlim(lims); ylim(lims);
    %cb vs cr greater than cluster
    a(end+1) = subplot(row,col,ii);ii=ii+1; plot(tg(:,cbi),tg(:,cmi),'ro'); title('rayleigh > 0.5 corr b vs m');
    hold on; plot([min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    xlim(lims); ylim(lims);
    %hists
    bn = -0.05:0.005:0.05; 
    a(end+1) = subplot(row,col,ii);ii=ii+1;hist(tl(:,cbi),bn);title(sprintf('%s\n%s','bef: hist corrs','~(r1&r2 > 0.5)'));xlim([-0.05 0.05]);
    a(end+1) = subplot(row,col,ii);ii=ii+1;hist(tg(:,cbi),bn);title(sprintf('%s\n%s','bef: hist corrs',' (r1&r2 > 0.5)'));xlim([-0.05 0.05]);
    a(end+1) = subplot(row,col,ii);ii=ii+1;hist(tl(:,cmi),bn);title(sprintf('%s\n%s','mid: hist corrs','~(r1&r2 > 0.5)'));xlim([-0.05 0.05]);
    a(end+1) = subplot(row,col,ii);ii=ii+1;hist(tg(:,cmi),bn);title(sprintf('%s\n%s','mid: hist corrs',' (r1&r2 > 0.5)'));xlim([-0.05 0.05]);
    %end
    
     %RAYLEIGH UIX
    f = figure(2); tp = uix.HBox( 'Parent', f, 'Spacing', 0);lp = uix.VBox('Parent', tp );rp = uix.VBox('Parent', tp);
    set( tp, 'Widths', [-1,-1] );% Adjust the main layout
    %vc= uicontainer('Parent', lp);% t = [t -1];
    row = 4; col = 2; ii =1; a = [];
    t = iud; lims = [min(min(t(:,cbi), min(t(:,cmi)))) max(max(t(:,cbi)), max(t(:,cmi)))] ;
    tg = t(t(:,r2mi)>=0.5 & t(:,r1mi)>=0.5,:); tl = t(t(:,r2mi)<0.5 | t(:,r1mi)<0.5,:);  
    %r1 vs r2
    a(end+1) = axes('Parent',uicontainer('parent',lp)); plot(a(end),t(:,r1mi),t(:,r2mi),'go'); 
    xlim(a(end),[0 1]); ylim(a(end),[0 1]); title(a(end),sprintf('decreasing gs pairs (%.2f>%.2f)\nrayleigh i1 vs i2 MID',gridThreshBef,gridThreshMid));
    %cb vs cr less than cluster
    a(end+1) = axes('Parent',uicontainer('parent',lp)); plot(a(end),tl(:,cbi),tl(:,cmi),'ro'); title(a(end),'rayleigh < 0.5 corr b vs m');
    hold(a(end),'on'); plot(a(end),[min(tl(:,cbi)) max(tl(:,cbi))], polyval( polyfit( tl(:,cbi),tl(:,cmi),1),[min(tl(:,cbi)) max(tl(:,cbi))] ) ); 
    xlim(a(end),lims); ylim(a(end),lims);
    %cb vs cr greater than cluster
    a(end+1) = axes('Parent',uicontainer('parent',lp)); plot(a(end),tg(:,cbi),tg(:,cmi),'ro'); title(a(end),'rayleigh > 0.5 corr b vs m');
    hold(a(end),'on'); plot(a(end),[min(tg(:,cbi)) max(tg(:,cbi))], polyval( polyfit( tg(:,cbi),tg(:,cmi),1),[min(tg(:,cbi)) max(tg(:,cbi))] ) );
    xlim(a(end),lims); ylim(a(end),lims);
    %hists
    bn = -0.05:0.005:0.05; %'Padding', 0, 'Spacing', 0 
    rp1 = uicontainer('Parent', rp);%title(a(end),sprintf('%s\n%s','bef: hist corrs','~(r1&r2 > 0.5)'));title(a(end),sprintf('%s\n%s','bef: hist corrs',' (r1&r2 > 0.5)'));
    a(end+1) = axes('Parent',uicontainer('parent',rp1,'position',[0,0  ,1,0.45])); hist(a(end),tg(:,cbi),bn);xlim(a(end),[-0.05 0.05]);
    alpha(a(end),0.9);  box(a(end),'on'); 
    a(end+1) = axes('Parent',uicontainer('parent',rp1,'position',[0,0.38,1,0.45])); hist(a(end),tl(:,cbi),bn);xlim(a(end),[-0.05 0.05]); set(a(end),'xticklabel',[]);%xlabel(a(end),[]); 
    alpha(a(end),0.9);box(a(end),'on'); 
    rp2 = uix.VBox('Parent', rp,'Padding', 0, 'Spacing', 0 );%title(a(end),sprintf('%s\n%s','mid: hist corrs','~(r1&r2 > 0.5)'));title(a(end),sprintf('%s\n%s','mid: hist corrs',' (r1&r2 > 0.5)'));
    a(end+1) = axes('Parent',uicontainer('parent',rp2)); hist(a(end),tl(:,cmi),bn);xlim(a(end),[-0.05 0.05]);set(a(end),'xticklabel',[]);%axis(a(end),'xticklabel',[]); 
    a(end+1) = axes('Parent',uicontainer('parent',rp2)); hist(a(end),tg(:,cmi),bn);xlim(a(end),[-0.05 0.05]);
    set( rp, 'heights', [-1,-1] );
    %end
    
    
        % MODULE
    criuBef = {}; criuMid = {};
    criuB =[]; criuM = [];
    for ri = 1:length(gids)
        gis = gids(ri)
        cis = cids{gis}
        g = groups{gis};
        g = g(cis);
        for j = 1:length(g)-1
            for k = j+1:length(g)
                c1 = g(j).before; c2 = g(k).before;
                c = shuffleTimeCorrelations (c1,c2, 1, 200, 0.006); 
                r1 = c1.rayleigh_score; r2 = c2.rayleigh_score;
                iu = Intersection_Over_Union(c1.module.x,c1.module.y, c2.module.x,c2.module.y);
                criuBef{g(j).ind,g(k).ind} = [c r1 r2 iu];
                criuB(end+1,:) = [g(j).ind g(k).ind c r1 r2 iu];
                c1 = g(j).midall; c2 = g(k).midall;
                c = shuffleTimeCorrelations (c1,c2, 1, 200, 0.006);
                r1 = c1.rayleigh_score; r2 = c2.rayleigh_score;
                iu = Intersection_Over_Union(c1.module.x,c1.module.y, c2.module.x,c2.module.y);
                criuMid{g(j).ind,g(k).ind} = [c r1 r2 iu];
                criuM(end+1,:) = [g(j).ind g(k).ind c r1 r2 iu];
            end
        end
    end

    figure; %hist iu
    subplot(1,3,1);hist(iud,20); title('IU pairs decreasing'); xlabel('intersection / union');
    subplot(1,3,2);hist(iun(iun>0),20); title('IU pairs non-decreasing');
    subplot(1,3,3);hist(iua(iua>0),20); title('IU pairs all'); 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SHUFFLING
    gcorrBefore = {}; gcorrMid = {};
    corrBefore = [];  corrMid = [];
    p.n = 250; p.movmean = 32;p.sigma = 0; p.bandf = -1; p.lag = 0;
    tic
    for ri = 1:length(gids)
        gis = gids(ri);
        cis = cids{gis};
        g = groups{gis};
        g = g(cis);
        tic
        for j = 1:length(g)-1
            for k = j+1:length(g)
                [ri j k] 
                gcorrBefore{k,j,ri} = shuffleTimeCorrelations (g(j).before, g(k).before,p);
                gcorrMid{k,j,ri}    = shuffleTimeCorrelations (g(j).midall, g(k).midall,p);
                corrBefore{g(j).ind,g(k).ind} = gcorrBefore{k,j,ri};
                   corrMid{g(j).ind,g(k).ind} = gcorrMid{k,j,ri};
            end
        end
        toc
    end
    disp('************************');
    toc
      
    

    ri = 11, ci = 12, c1 = g(1).before, c2 = g(2).before, 
    [r, c] = size(corrBefore);
    pvalsBefMid = cell(r,c); pbs = []; pms = []; ijcbmpbm = []; 
    for ri = 1:r
        for ci = 1:c
            if ~isempty(corrBefore{ri,ci})
                t = corrBefore{ri,ci}; cb = t(1);
                [s,I] = sort(abs(t),'descend');
                pb = round(find(I==1)/length(I),2); %index of non shuffled
                t = corrMid{ri,ci}; cm = t(1);
                [s,I] = sort(abs(t),'descend');
                pm = round(find(I==1)/length(I),2);
                pvalsBefMid{ri,ci} = [pb pm]; pbs(end+1) = pb; pms(end+1) = pm;
                ijcbmpbm(end+1,:) = [ri ci cb cm pb pm ];
            end
        end
    end
    figure(1);
    subplot(2,1,1);hist(pbs*100,100); title('p shuffle before 250 update');
    subplot(2,1,2);hist(pms*100,100); title('p shuffle mid 250 update');
    
    i1i = 1; i2i = 2; cbi = 3; cmi = 4; pbi = 5; pmi = 6; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FIGURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f = figure;%tp = uix.HBox( 'Parent', f, 'Spacing', 0); lp = uix.VBox('Parent', tp );rp = uix.VBox('Parent', tp);
 
    %2 HISTS C B M
    figure; 
    bn = [-.05:0.0025:.05]; %LEARN DIFFERENCES BETWEEN HIST AND HISTOGRAM
    subplot(1,2,1);histogram(ijcbmpbm(:,cbi),bn); hold on
    %title(sprintf('%s\n%s','Time Correlation Between Simultaneously Recorded Grid Cells','Pre Muscimol'),'fontsize',18);
    title(sprintf('%s','Pre Muscimol'),'fontsize',18); 
    histogram(ijcbmpbm(ijcbmpbm(:,pbi)<=0.01,cbi),bn,'FaceColor','r','EdgeColor','r');
    xlabel('Correlation Value','fontsize',16);
    ylim([0,35]);
    subplot(1,2,2);histogram(ijcbmpbm(:,cmi),bn); title('During Muscimol','fontsize',18); hold on
    histogram(ijcbmpbm(ijcbmpbm(:,pmi)<=0.01,cmi),bn,'FaceColor','r','EdgeColor','r');
    xlabel('Correlation Value','fontsize',16);
    legend({sprintf('%s\n%s','corrs for pairs passing','gridscore threshold'),...
        sprintf('%s\n%s','P-val above 0.01','of Shuffling')},'fontsize',16);
    
    %SPARSE
    
    %SCATTER CORR CORR P-Value
    h= figure(2);clf(h);set(h,'color','w');  afs = 10; mfs = 10; tfs = 12; lfs = 3;
    plot(ijcbmpbm(:,cbi),ijcbmpbm(:,cmi),'.','markersize',mfs,'linewidth',lfs); hold on;
    pd.shufbx=ijcbmpbm(ijcbmpbm(:,pmi)<=0.01,cbi);pd.shufmy=ijcbmpbm(ijcbmpbm(:,pmi)<=0.01,cmi);
    plot(     ijcbmpbm(ijcbmpbm(:,pmi)<=0.01,cbi),          ijcbmpbm(ijcbmpbm(:,pmi)<=0.01,cmi),...
         'mo','markersize',mfs,'linewidth',lfs);
     %plot([min(ijcbmpbm(:,cbi)) max(ijcbmpbm(:,cbi))], polyval( polyfit( ijcbmpbm(:,cbi),ijcbmpbm(:,cmi),1),[min(ijcbmpbm(:,cbi)) max(ijcbmpbm(:,cbi))]));
     f1=fit(ijcbmpbm(:,cbi),ijcbmpbm(:,cmi),'poly1');
     pscat = plot(f1,'b'); pscat.LineWidth = lfs;
     legend({'corr pre vs dur',...
         sprintf('%.2fx + %.2f',round(f1.p1,2),round(f1.p2,2)),...
         'significance (dur) < %1'...
        },'fontsize',afs,'location','best');
    %ylim([-.02 .03]); xlim([-.02 .03]); 
    xlabel('correlation before','fontsize',afs); 
    ylabel('correlation during','fontsize',afs)
    title('Time Correlation Pre vs During by Shuffling P-value','fontsize',tfs);
    %set(gca,'XTickLabel',[-0.02:0.01:0.03],'fontsize',14)
    %set(gca,'YTickLabel',[-0.02:0.01:0.03],'fontsize',14)
    axis square; axis equal;
   
    figure
    f1=fit(ijcbmpbm(:,cbi),ijcbmpbm(:,cmi),'poly1');%fit(cdate,pop,'poly1')
    plot(f1);
    legend(sprintf('%.2fx+%.2f',round(f1.p1,2),round(f1.p2,2)));
    
    disp('');
    
    
    params.lag_max_secs = 2; params.sigma = 1; params.number_degree_bins = 1;
    str.tt = '';str.xt = '';str.yt = '';
    figure('Position', [50,50, 1250, 10000]);     %set(gcf, 'PaperSize', [20 20])    % Same, but for PDF output
    
    figure(1); r = 9; c = 8; clf; set(gcf, 'PaperPosition', [0 0 7 11])    % can be bigger than screen 
    l = r%length(iud)
    for i = 1:length(iud)
        %figure
        a = cells{iud(i,1)}; b = cells{iud(i,2)};     
        if mod(i-1,r) == 0; clf; end
        %BEFORE
        c1 = a.before;  c2 = b.before; %str.yt = 'PRE'; %ax = subplot(r,c,ii+1); 
        ii = mod(i-1,r)*c; p = 'position';
        ax = axes(p,apos(c,r*c,ii+1)); str.tt = sprintf('PRE: c%i',a.ind);plotRM(c1,str,ax);axis(ax,'off');
        ax = axes(p,apos(c,r*c,ii+2)); str.tt = sprintf('C%i',b.ind);plotRM(c2,str,ax);axis(ax,'off');
        ax = axes(p,apos(c,r*c,ii+3)); imagesc(ax, xcorr2(c2.rm,c1.rm));set(ax,'YDir','normal'); title('space corr');
        axis(ax,'square');colormap(ax,'jet'); axis(ax,'off'); str.tt ='time corr';
        ax = axes(p,apos(c,r*c,ii+4)); plotTimeCorr(c1,c2,params,str,ax); ax.YAxis.Visible = 'off';
        %DUR
        c1 = a.midall;  c2 = b.midall; %str.yt = 'DUR';
        ax = axes(p,apos(c,r*c,ii+5)); str.tt = sprintf('DUR: c%i',a.ind);plotRM(c1,str,ax);axis(ax,'off');str.yt = '';
        ax = axes(p,apos(c,r*c,ii+6)); str.tt = sprintf('C%i',b.ind);plotRM(c2,str,ax);axis(ax,'off');
        ax = axes(p,apos(c,r*c,ii+7)); imagesc(ax, xcorr2(c2.rm,c1.rm));set(ax,'YDir','normal');title('space corr');
        axis(ax,'square');colormap(ax,'jet'); axis(ax,'off');  str.tt ='time corr';
        ax = axes(p,apos(c,r*c,ii+8)); plotTimeCorr(c1,c2,params,str,ax); ax.YAxis.Visible = 'off';
        if mod(i-1,r) == r-1; print(gcf, sprintf('a%d.png',ceil(i/r)), '-dpng', '-r300' ); end
    end
    print(gcf, sprintf('a%d.png',ceil(i/r)), '-dpng', '-r300' ); i
    %}
    
    %SINGLE RAYLEIGH
    cids = unique([tl(:,1) tl(:,2); tg(:,1) tg(:,2)]);
    for i = 1:length(cids);
        j = cids(i);
        rs(i,:) = ...
            [j cells{j}.before.rayleigh_score cells{j}.midall.rayleigh_score]
    end
    
    clf
    plot(rs(:,2),rs(:,3),'b.','markersize',mfs);
    hold on; plot(rs(rs(:,3)>0.5,2),rs(rs(:,3)>0.5,3),'r.','markersize',mfs);
    %title('rayleigh score by cell');
    legend({'Cluster1';'Cluster2'},'fontsize',afs,'Location','best');
    set(gca,'XTick',[0,0.5,1]);
    set(gca,'YTick',[0,0.5,1]);
    set(gca,'XTickLabel',[0,0.5,1],'fontsize',14)
    set(gca,'yTickLabel',[0,0.5,1],'fontsize',14)
    xlabel('before','fontsize',afs);
    ylabel('during','fontsize',afs);
    %tg = t(t(:,r2mi)>=0.5 & t(:,r1mi)>=0.5,:);
end





