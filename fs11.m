%f3s6    3s6
function fs11(cellsn,pairs)       
        cang=[];
        load('./data/fs11_cang');
%         for j=1:len(pairs)
%             j
%             c = cellsn(pairs(j,1)); d = cellsn(pairs(j,2));
%             [~, a]=gridscore2(xcorr2g(c.before.rm,d.before.rm),2); %xcorr2g xcorr2n
%             [~, b]=gridscore2(xcorr2g(c.midall.rm,d.midall.rm),2);
%             [~, e]=gridscore2(xcorr2g(c.after.rm, d.after.rm) ,2);
%             %[wrapTo180(a) wrapTo180(b) wrapTo180(e)]
%             cang(j,:) = [wrapTo180(a) wrapTo180(b) wrapTo180(e)]; %[deg2rad(a) deg2rad(b) deg2rad(e)];
%         end    
        all=wrapTo180(cang(:,2)-cang(:,1));
        all2=wrapTo180(cang(:,3)-cang(:,1));
        
        figure(1011); clf; set(gcf,'color','w','position',[10 10 900 800]); 
        %hists
        subplot(222)
        bar(sort(abs(all))); b= gca(); %b.YLim=[0 180];
        b.XLabel.String = 'pair'; b.YLabel.String='angle';
        b.Title.String = 'difference in angle pre vs dur';
        subplot(224)
        bar(sort(abs(all2))); b= gca(); %b.YLim=[0 180];
        b.XLabel.String = 'pair'; b.YLabel.String='angle';
        b.Title.String = 'difference in angle pre vs post';
        subplot(221); hold on;  axis square;
        %scatter
        d1=cang(:,1);d2=cang(:,2); plot(d1,d2,'.','markersize',10);
        b= gca(); b.XLabel.String = 'pre'; b.YLabel.String='dur';
        b.Title.String = 'first field angle pre vs dur';
        [r,p]=circ_corrcc(deg2rad(d1),deg2rad(d2));
        %text(0.03,0.98,sprintf('r=%.3f p=%.3f',round(r,3), round(p,3)),...
        text(0.03,0.98,sprintf('r=%.3f %s',round(r,3), pstr(p,3)),...
            'Units','normalized','color','r');
        subplot(223); hold on;  axis square;
        d1=cang(:,1);d2=cang(:,3); plot(d1,d2,'.','markersize',10);
        b= gca(); b.XLabel.String = 'pre'; b.YLabel.String='post';
        b.Title.String = 'first field angle pre vs post';
        [r,p]=circ_corrcc(deg2rad(d1),deg2rad(d2));
        %text(0.03,0.98,sprintf('r=%.3f p=%.3f',round(r,3), round(p,3)),...
        text(0.03,0.98,sprintf('r=%.3f %s',round(r,3), pstr(p,3)),...
            'Units','normalized','color','r');       
        suptitle('Cell pair spatial cross-correlation: angle of first firing field');
 
    end 