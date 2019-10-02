%{
load('C:\\Noam\\Data\\muscimol\\noam\\cells_45min_d_patchtraj_rayleigh'); %personal

load('C:\\Noam\\Data\\muscimol\\cells15nan')
load('.\data\pairs');
cels = unique([pairs(:,1),pairs(:,2)])';
%}
function Revision(cellsn, pairs)




dbstop if error

cels = unique(pairs(:))';
clls = [cellsn(cels)]; %aclls=cellsn;
s = {'before';'midall';'after'};tss = {'pre','dur','post'};
vs='rayleigh_score';va='rayleigh_angle';

 reprocess
 %fieldAngle
 %hdangle;
 %drift;


    function reprocess()
        a=[];b=a;c=a;d=a;e=a;f=a;g=a;h=a;m=a;n=a;q=a;t=a;w=a;x=a;y=a;z=a; i=a;     
        
        
        b=[cellsn.before]; l =     arrayfun(@(z) len(z.st)>=100,b); 
        m=[cellsn.midall]; l = l & arrayfun(@(z) len(z.st)>=100,m); 
        b=[b.gridscore]; m=[m.gridscore];
        g= b>0 & m<0; clear b; clear m;
        t=cellsn(g&l); clear g; clear l;
        t=[t.ind];

        %t=unique(pairs(:));
        figure(1);
        for i=t
            clf;
            c=cellsn(i).before;
            subplot(241)
            imgsc(c.rm,1); title('before');
            subplot(243)
            imgsc(c.ac3,2); title(['ac3 ' n2(c.gs3)]);
            %c=cellsn(i).midall;
            subplot(242)
            rm=interpnanrm(c.rm);
            %imgsc(normxcorr2ni(rm),2); title(['interp' n2(gridscore2(normxcorr2ni(rm),2))]);
            %imgsc(c.ac2,2); title(['ac2 ' n2(c.gs2)]);
            subplot(244)
            imgsc(c.ac,2); title(['ac ' n2(c.gridscore)]);
            suptitle(n2(i));
                        c=cellsn(i).midall;
            subplot(245)
            imgsc(c.rm,1); title('midall');
            subplot(246)
            %imgsc(c.ac3,2); title(['ac3 ' n2(c.gs3)]);
            %c=cellsn(i).midall;
            subplot(247)
            rm=interpnanrm(c.rm);
            imgsc(c.ac3,2); title(['ac3 ' n2(c.gs3)]);
            %imgsc(normxcorr2ni(rm),2); title(['interp' n2(gridscore2(normxcorr2ni(rm),2))]);
            %imgsc(c.ac2,2); title(['ac2 ' n2(c.gs2)]);
            subplot(248)
            imgsc(c.ac,2); title(['ac ' n2(c.gridscore)]);
            suptitle(n2(cellsn(i).ind));
        end
        
        assert(false);
        
%        t = aclls; tl=arrayfun(@(z) len(z.st),t); 
        
        %nancllg=cellsn;
        s = {'before';'midall';'after'};tss = {'pre','dur','post'};
        %changing ratemap to nans noam method
        for ss=[s(:)']
            for i = 1:len(cellsn)
                [n2(i) ss]
                c=cellsn(i).(ss{:});z=c; 
%                 [c.rm, c.max_r] = createRateMapNan(c,50);
%                 %cellsn(i).(ss{:}).max_r=max_r;
%                 c.ac = xcorr2g(c.rm,c.rm); <IMPORTANT DO NOT DELETE

                  ac3=normxcorr2ni(c.rm);
                  cellsn(i).(ss{:}).ac3=ac3;
                  cellsn(i).(ss{:}).gs3=gridscore2(ac3,2);
                  
                  ac2=xcorr2n(c.rm);
                  cellsn(i).(ss{:}).ac2=ac2;
                  cellsn(i).(ss{:}).gs2=gridscore2(ac2,2);
%                 c.gridscore=gridscore2(c.ac, 2);
%                 [c.gridscore z.gridscore]
%                 c.module = [];
%                 c.module = Find_Module(imgaussfilt(c.ac, 3,'FilterDomain','spatial'));
%                 nancllg(i).(ss{:})=c; 
            end
        end
        
        start the simulation on the lab computer
        
        d=[aclls.  ('before')]; d=[d.gridscore]; e=[aclls.  ('midall')]; e=[e.gridscore];
        f=[nanclls.('before')]; f=[f.gridscore]; g=[nanclls.('midall')]; g=[g.gridscore];
        j=[nancllg.('before')]; j=[j.gridscore]; k=[nancllg.('midall')]; k=[k.gridscore]; 
        clf; plot(d,f,'.');axis equal;
        clf; plot(d,j,'.');axis equal;
        clf; plot(e,k,'.');axis equal;
        sum(d>0.3 & e<0.25) 
        sum(f>0.3 & g<0.25)
        sum(j>0.3 & k<0.25)
        
        a=cellsn(111).before.rm;
        %a(isnan(a))=0; a(isnan(a))=nanmean(a(:)); 
        b=a;
        c=xcorr2g(a,b);
        imgsc(c,0.1,2);
        imgsc(abs(c),0.1,2);
        imgsc(c,0.1,2);
        
        
        
        %making midall 15 min or later
        for e = 1:len(aclls)
            c=aclls(e).('midall'); z=t(e);
            e
            a=c.pt>=15*60; d=[]; [c.pt(1) len(c.pt) sum(a)]
            d.pt=c.pt(a);d.px=c.px(a);d.py=c.py(a);d.hd=c.hd(a);
            a=c.st>=15*60;
            d.st=c.st(a);d.sx=c.sx(a);d.sy=c.sy(a);
            %rm
            c=d;
            [d.rm d.max_r]=createRateMap(c, 50);
            %figure(2);subplot(211); imgsc(z.rm);subplot(212); imgsc(d.rm);
            [d.rayleigh_score, d.rayleigh_angle ,d.hd, d.shd]=...
                rayleigh_score(c.pt,c.st,c.hd);
            [rad2deg(z.rayleigh_angle) rad2deg(d.rayleigh_angle)]
            
            d.ac = xcorr2(d.rm); 
            d.gridscore = gridscore2(d.ac, 2);
            [d.gridscore z.gridscore]
            d.module = Find_Module(imgaussfilt(d.ac, 3,'FilterDomain','spatial'));
            d.exists = false;
            if d.max_r ~= 0 
                d.exists = true;
            end
            aclls(e).('midall')=d;
        end        
        
        z = [aclls.('midall')];
        y=arrayfun(@(y) len(y.st), z);     
        a=[z.rayleigh_score]; b=[t.rayleigh_score]; sum(a>0.3&y>=100)
        sum(b>0.3&tl>=100)
        c=[z(cels).gridscore];  [~, x]=max(c); 
        m=t(cels(x)); n=z(cels(x));
        figure(2);
        subplot(211);imgsc(m.rm,3); title(m.gridscore); sum(m.rm,'all')
        subplot(212);imgsc(n.rm,3); title(n.gridscore); sum(n.rm,'all')
        gridscore2(xcorr2(n.rm), 2)
        
    end %reprocess

    


    function drift()
        a=[];b=a;c=a;d=a;e=a;f=a;g=a;h=a;m=a;n=a;q=a;t=a;w=a;x=a;y=a;z=a;
        
         
        nbins = 100; twins = [1 2 3 5 10] ; %window in secs
        s = {'before', 'midall','after'}; ts = {'pre','dur','post'};
        
        t=[aclls.(s{1})];d=cellfun(@(x) max(diff(x)),{t.pt},'uni',1);
        t=[aclls.(s{2})];e=cellfun(@(x) max(diff(x)),{t.pt},'uni',1); off=max([d e]);%0.04
        
        tic
        load('.\\data\\dxdyrate','dxywinrd','dxywinrn','gridwin');%,'autowin');
        toc
        
  % vector sum of angle: add all sins and all coss of each angle, arctan of all; or mean sins / mean cos
        
        
%         figure(1); colormap jet;
%         subplot(211); imagesc(dxywinrd{93,5}.(s{1})); title('dmet');  axis square; 
%         subplot(212); imagesc(dxywinrn{93,5}.(s{1})); title('nmet');  axis square;
         
        e=dxywinrd{16,10}.(s{1});  g=Cross_Correlation(e,e); 
        %f=ones(size(e)); y=xcorr2(f,f); m=max(y(:));  y=ones(size(y))*m; % numerator        
        %f(isnan(e))=0; z=xcorr2(f,f); z=y./z; z(z>100)=0; x=h.*z;  
        x=xcorr2n(e); h=e; h(isnan(h))=0; h=normxcorr2(h,h); % x=e; x(isnan(x))=0; x=normxcorr2(x,x);
        a=2; d=99; d = 100-d:100+d; 
        figure(1); clf; colormap jet;
        subplot(221); imgsc(e,a); title('input');
        %subplot(332); imgsc(f); title('non-nans');
        subplot(223); imgsc(g(d,d),a); title(['gmet ' n2(round(gridscore2(g,a),2)) ]); 
        subplot(222); imgsc(h(d,d),a); title(['normxcorr2 ' n2(round(gridscore2(h,a),2)) ]); 
        subplot(224); imgsc(x(d,d),a); title(['nmet ' n2(round(gridscore2(x,a),2)) ]);
        
      
       
        %gridwin={};
        %histogram scatters
        dxyr=dxywinrd; z='dmet';%dxyr=dxywinrn; z='nmet';
        m=1;w=2; mm=['corr ' ts{m} ' ' ts{w}]; g=[];%compare against
        figure(3000+w); clf; pr = [200 10 900 600]; set(gcf,'position',pr);
        for h =1:len(twins)
            twin=twins(h);
            g={}; %scatter
            for w = 2:3 %scatter
                a=[];
                for i = 1:len(dxyr)
                %tic; [twin i]
                %corr sesh vs sesh for dxdy
                t=corr(dxyr{i,twin}.(s{m})(:),dxyr{i,twin}.(s{w})(:),'rows','complete');
                %gridscore sesh sesh dxdy
                %to make e=dxyr{i,twin}.(s{1}); f=dxyr{i,twin}.(s{2}); d=dxyr{i,twin}.(s{3}); 
                %   e=xcorr2nan(e,e); f=xcorr2nan(f,f);d=xcorr2nan(d,d);autowin{i,twin}.before=e; 
                %   autowin{i,twin}.midall=f; autowin{i,twin}.after=d;
                %t=[gridscorenan(e,5) gridscorenan(f,5)]
                a(end+1,:)=t; %toc
                end
            %g(:,h)=a; %g{h}=a(~isnan(a));%
            g{w}=a; %scatter
            %gridwin{twin}=a;
            end
           subplot(2,3,h);
            %scatter
            y=toCol(g{2});x=toCol(g{3});t=~isnan(y)&~isnan(x);x=x(t);y=y(t); scatter(x,y); hold on;
            f = fit(x, y,'poly1');plot(x,f(x),'-'); axis('tight'); axis square; [a b]=ccof(x,y);
            text(0.1,0.9,sprintf('a=%.2f r=%.2f p=%.3f',round(f.p1,2), round(a,2), round(b,3)),'Units','normalized');
            title([n2(twin) 's']); xlabel('pre vs post'); ylabel('pre vs dur'); 
            e=slim(gca);xlim(e);ylim(e);
        end
        %scatter
        t=['correlation of dxdy by time window ' z];
        suptitle(t);
        saveas(gcf,['./figs/scatter ' t '.png']);     




a=[1:3;1:3]


        z='dmet';
        m=1;w=2;t=['gridscore nan of dxdy by time window sigma-5 ' z];%compare against
        figure(5000); clf; pr = [200 10 800 600]; set(gcf,'position',pr);
        for h =1:len(twins)
            twin=twins(h);           
            subplot(2,3,h);
            g=gridwin{twin};
            %scatter
            y=g(:,2);x=g(:,1);scatter(x,y); hold on;
            f = fit(x, y,'poly1');plot(x,f(x),'-'); axis('tight'); axis square;
            text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f.p1,3), round(ccof(x,y),3)),'Units','normalized');
            title([n2(twin) 's']); xlabel(ts{m}); ylabel(ts{w}); 
            e=slim(gca);xlim([-2 2]);ylim([-2 2]);
        end
        suptitle(t);
        saveas(gcf,['./figs/' t '.png']);%plot(sort(a));  


   a;     
        %DMETH
        %{
        %dxywinrd={};
        for h = 1:len(twins)
            twin = twins(h);
            for j=1:len(pairs)
                ['a ' n2(twin) ' ' n2(j)]
                tic
                for ss=[s(:)']
                    a = aclls(pairs(j,1)).(ss{:}); b = aclls(pairs(j,2)).(ss{:});
                    rmsa = nan(nbins); rmta=rmsa; tssi = 1;tppi = 1; %norm=zeros(nbins);
                    for i = 1:len(a.st)
                        sr = driftwin(twin,a.st(i),a.sx(i),a.sy(i),b.st,b.sx,b.sy,tssi,0);
                        pr = driftwin(twin,a.st(i),a.sx(i),a.sy(i),b.pt,b.px,b.py,tppi,off);
                        %z=~ismember(sr.a(:),pr.a(:)); assert(sum(z)==0);
                        rms = driftmap(sr.a,max(a.px),max(a.py)); tssi = sr.i;
                        rmt = driftmap(pr.a,max(a.px),max(a.py)); tppi = pr.i;
                        %there should not be bins with spike and no time.
                        rmt(~isnan(rms) & isnan(rmt)) =rms(~isnan(rms) & isnan(rmt)); %ADD BACK
                        %should not be bins with more spikes than time?
                        rmt(rms>rmt)=rms(rms>rmt);
                        t = nansum(cat(nbins,rmsa,rms),nbins); 
                        t(isnan(rmsa) & isnan(rms))=nan; %only maintain nan if both values are not nan.
                        rmsa=t;
                        t = nansum(cat(nbins,rmta,rmt),nbins); %sets all nans to 0's 
                        t(isnan(rmta) & isnan(rmt))=nan; %only keep nan if both values are nan.
                        rmta=t;
                    end
                    
                    %assert(isequal(norm==0,isnan(dxyr)));
                    dxyr=rmsa./rmta;
                    %assert( max(dxyr(:) ) <= 1 );
                    dxywinrd{j,twin}.(ss{:}) = dxyr;
                end
                toc;
            end
        end
        %}
        
        
        
        %NMETH
        %{
        %dxywinrn={};
        for h = 1:len(twins)
            twin = twins(h);
            for j=1:len(pairs)
                ['n ' n2(twin) ' ' n2(j)]
                tic
                for ss=[s(:)']
                    a = aclls(pairs(j,1)).(ss{:}); b = aclls(pairs(j,2)).(ss{:});
                    dxyr = nan(nbins);tssi = 1;tppi = 1; norm=zeros(nbins);
                    for i = 1:len(a.st)
                        sr = driftwin(twin,a.st(i),a.sx(i),a.sy(i),b.st,b.sx,b.sy,tssi,0);
                        pr = driftwin(twin,a.st(i),a.sx(i),a.sy(i),b.pt,b.px,b.py,tppi,off);
                        %z=~ismember(sr.a(:),pr.a(:)); assert(sum(z)==0);
                        rms = driftmap(sr.a,max(a.px),max(a.py)); tssi = sr.i;
                        rmt = driftmap(pr.a,max(a.px),max(a.py)); tppi = pr.i;
                        assert(sum(~isnan(rms(:)) & isnan(rmt(:)))==0);
                        rmt(~isnan(rms) & isnan(rmt)) =rms(~isnan(rms) & isnan(rmt));
                        %assert(sum(rms(:)>rmt(:))==0);
                        %if sum(rms(:)>rmt(:)) > 0
                        %['s > t ' n2(sum(rms(:)>rmt(:)))]
                        %end
                        rmt(rms>rmt)=rms(rms>rmt); %How does this happen?
                        %rmt(rmt==0) = inf; BAD Noam
%                       why spike bin where there is no movement??
                        t = rms./rmt; %t(t==inf)=nan; %now impossible
                        tt=dxyr; %0/0 = nan, should I be getting any infs???
                        dxyr = nansum(cat(nbins,tt,t),nbins); %dxyr+t
                        dxyr(isnan(t) & isnan(tt))=nan; %only maintain nan if both values are not nan.
                        norm=norm + ~isnan(t); %add a value each time value is not nan to normalize.
%                         figure(1);
%                         subplot(221); imagesc(rms);set(gca,'ydir','normal');
%                         subplot(222); imagesc(rmt);set(gca,'ydir','normal');
%                         subplot(223); imagesc(dxyr);set(gca,'ydir','normal');
%                         subplot(224); imagesc(norm);set(gca,'ydir','normal');
                    end
                    %assert(isequal(norm==0,isnan(dxyr)));
                    dxyr=dxyr./norm;
                    %assert( max(dxyr(:) ) <= 1 );
                    dxywinrn{j,twin}.(ss{:}) = dxyr;
                end
                toc;
            end
        end
        %}

        
        
        % % % % graphics imagesc 
        d = 16; d = 50-d:51+d; %b=1;
        pr = [200 10 1800 900]; s = {'before','midall'};  t=2; e = 2;% f=5;f=ceil(sqrt(f));%f-#pairs %t 
        for twin = 1%twins
            figure(1000+10*twin+t); clf; set(gcf,'position',pr);
            for i = 1:20%(pairs)
                subplot(4,e*5,e*i-e+1);x = 1; %PRE 1 DUR 2
                a = dxywinrd{i,twin}.(s{x}); a=a(d,d); %b = gridscore2(a,t);
                imagesc(imgaussfilt(a, t,'FilterDomain','spatial'));
                %imagesc(xcorr2(imgaussfilt(a, t,'FilterDomain','spatial')));
                %a(isnan(a))=-1;%imagesc(a);
                colormap jet; axis off; axis square;title(sprintf('p%d:%s',i,ts{x}));
                
                subplot(4,e*5,e*i); x = 2;%PRE 1 DUR 2
                a = dxywinrd{i,twin}.(s{x}); a=a(d,d); % b = gridscore2((a),t);
                imagesc(imgaussfilt(a, t,'FilterDomain','spatial'));
                %imagesc(imgaussfilt(xcorr2(a), t,'FilterDomain','spatial'));
                colormap jet; axis off; axis square;title(sprintf('p%d:%s',i,ts{x}),'color','red');
            end
            %suptitle(sprintf('%s %ds: dxdy-spike/dxdy-pos and gridscore by pair sigma=%d',s{sr},twin,t));
            suptitle(sprintf('%ds: dxdy-spike/dxdy-pos and gridscore by pair sigma=%d',twin,t));
            saveas(gcf, sprintf('./figs/drift_diff_%ds_sig%d_zoom_dmet.png',twin,t) );
        end
        
        % % % % Time window lapse
        t=3;  d = 25; d = 51-d:51+d; b=1:20; %b = b+len(b);
        figure(1880+t); clf; pr = [100 10 1900 1000]; set(gcf,'position',pr);
        for e = 1:2
            for i = 1:len(twins)
                for j = 1:len(b)/2
                    c = (e-1)*len(b)/2+j; %which pair
                    subplot(len(twins)*2,len(b),  ((e-1)*len(twins)+i-1)*len(b)+2*j-1);sr = 1; cla;
                    a = dxywinr{b(c),twins(i)}.(s{sr}); a=a(d,d);
                    hold on; title(sprintf('P%d%s %ds',b(c),s{sr}(1:3),twins(i)));
                    if i>1 %subtract previous
                        f = dxywinr{b(c),twins(i-1)}.(s{sr}); f=f(d,d);
                        a=a-f;
                        title(sprintf('P%d%s %ds:%ds',b(c),s{sr}(1:3),twins(i)-twins(i-1),twins(i)) );
                    end
                    imagesc(imgaussfilt(a, t,'FilterDomain','spatial'));
                    colormap jet; axis off; axis square;
                    
                    subplot(len(twins)*2,len(b),  ((e-1)*len(twins)+i-1)*len(b)+2*j );sr = 2; cla
                    a = dxywinr{b(c),twins(i)}.(s{sr}); a=a(d,d);
                    hold on; title(sprintf('p%d%s %ds',b(c),s{sr}(1:3),twins(i)),'color','red');
                    if i>1 %subtract previous
                        f = dxywinr{b(c),twins(i-1)}.(s{sr}); f=f(d,d);
                        a=a-f;
                        title(sprintf('p%d%s %ds:%ds',b(c),s{sr}(1:3),twins(i)-twins(i-1),twins(i))...
                            ,'color','red');
                    end
                    imagesc(imgaussfilt(a, t,'FilterDomain','spatial'));
                    colormap jet; axis off; axis square;
                end
            end
        end
        suptitle(sprintf('row=window[1 3 5 10s] col=pair[pre dur] sigma=%d',t));
        
        % % % % gridscore fig
        figure(58); clf; hold on; t = 1;%sigma
        for h = 1:len(twins)
            %calculate all gridscores pre and dur
            twin = twins(h);
            a = [dxywinr{:,twin}]; a = {a.('before')};
            a = cellfun(@(x) gridscore2(xcorr2(x),t),a);%,'uni',0); %do I need this?
            b = [dxywinr{:,twin}]; b = {b.('midall')};
            b = cellfun(@(x) gridscore2(xcorr2(x),t),b);%,'uni',0);
            subplot(2,2,h)
            plot(a,b,'.');
            title(sprintf('drift %ds',twin)); axis square; axis equal;
            xlabel('gridscore pre'); ylabel('gridscore dur'); mm=slim(gca);
            xlim(mm);ylim(mm)
        end
        suptitle(sprintf('gridscore of dxdy-spike/dxdy-pos by drift window sigma=%d del %d:%d',t,del(1),del(end)));
        
        
        
      
        
        %%{ 

       
        dxywinro=dxywinr;
        rmt=dxywinr{11,10}.before;imagesc(imgaussfilt(rmt, 2,'FilterDomain','spatial'));
        
  
        %drift window
        %has rounding error between discrete x y vals and time x y
        function r = driftwin(twin,t,tx,ty,st,sx,sy,tsi,off)
            twin=twin+off;
            dxya = []; tri = tsi;
            emd=length(st); tmin=t-twin; tmax=t+twin;
            while tsi <= emd
                if   st(tsi) < tmin  % st(tsi) + twin < t
                    tri=tsi;%?????
                elseif  tmax < st(tsi) %t < st(tsi)-twin
                    tsi = emd;
                else
                    %if b.st(k)-twin <= t && t <= b.st(k)+twin
                    dx =  sx(tsi) - tx; %no need to round.
                    dy =  sy(tsi) - ty; %no need to round.
                    dxya(end+1,:) = [dx dy];
                end
                tsi=tsi+1;
            end
            r.a=dxya; r.i = tri;
            %assert (tri~=-1)
        end
        
        []
        %rate map
        function rm = driftmap(dxya,mpx,mpy) %a contains all diffs of x and y;
            rm = zeros(nbins);
            if length(dxya)>2
                dx = dxya(:,1);dy = dxya(:,2); %px = px+max(px); py = py+max(py);
                binx=-mpx:2*mpx/nbins:mpx;
                biny=-mpy:2*mpy/nbins:mpy;
                pxi = discretize(dx, binx);
                pyi = discretize(dy, biny);
                rm = rm + accumarray([pyi, pxi], 1, [nbins nbins]);
            else
                'empty drift map';
            end
            %rm(rm==0)=nan;
        end
     %}     
        
    end %drift()

    function hdangle()
        a=[],b=[],c=[],d=[],e=[],f=[],m=[],n=[],x=[],y=[],z=[];
        % ANGLE VILLE %%%% ^todo: alpha by score
        %aclls = [     cells{:}]; 
        clls = cellsn(cels);
        
        t = [clls.midall]; t = [t.rayleigh_score]; hdclls = clls(t>0.4);
        load('.\\data\\roomdose','room');
        %load('.\\data\\aclls15min','aclls');
        s = {'before';'midall';'after'}; vs='rayleigh_score';va='rayleigh_angle'; tss = {'pre','dur','post'};
        roomo=room;
        
        
        % RANGLE FOR ALL CELLS > RSCORE
        figure(702);clf;
        rst = 0.3;
        subplot(131); x = 1; t = [];fp(x,rst); len(t);
        subplot(132); x = 2; t = [];fp(x,rst); len(t);
        subplot(133); x = 3; t = [];fp(x,rst); len(t);
        suptitle(sprintf('r-angle for cells with r-score > %0.1f',rst));
        

        % vector sum of angle: add all sins and all coss of each angle, arctan of all; or mean sins / mean cos
        x=2; t = [cellsn.(s{x})]; tl=arrayfun(@(z) len(z.st),t);   
        ta = [t.(va)]; ts = [t.(vs)];  tis = ts>rst&tl>=100; ta = ta(tis);
        rad2deg(atan(sum(cos(ta))/sum(sin(ta))))
        
        
        function fp(x,rst)
            cla; xlabel('r angle');
            t = [cellsn.(s{x})]; tl=arrayfun(@(z) len(z.st),t); tl=tl>=100;   
            ta = rad2deg([t.(va)]); ts = [t.(vs)];  tis = ts>rst; t = ta(tis & tl);
            histogram(t,-180:30:180); xlim([-180 180]); ylim([0 40]); axis square;
            title(tss{x});
            text(0.1,0.95,sprintf('n=%d a=%.2f +- %.1f',len(t),mean(t),std(t)),'Units','normalized');
            %axis square;
        end        
        
        
        
        %room figure
        figure(2900);rst=0.3;st='all';
        %         figure(2901);clf;rst=0.4;st='grid2HD';
        clf; pr = [200 10 800 800]; set(gcf,'position',pr); f=[clls.(s{2})];
        for x=1:3%session
            t = [cellsn.(s{x})]; ts = [t.(vs)]'; tl=arrayfun(@(z) len(z.st),t)';room=roomo;
            %         t =  [clls.(s{x})]; ts = [f.(vs)]'; room=roomo; room=roomo(cels); %  [ HD ]
            z=unique(room);d= arrayfun(@(x) [t(ts>rst & room==x & tl>100).(va)],z,'uni',0);%z=unique(room(ts>rst));
            for y=1:len(z) %z to d
                subplot(3,len(z),(x-1)*len(z)+y);cla; %z to d
                a=rad2deg(d{y});histogram(a,-180:20:180);xlim([-180 180]); %ylim([0 len(a)]);
                xlabel('r angle'); axis square;
                title(sprintf('%s rm%d n=%da=%.0f+-%.0f',tss{x},z(y),len(a),mean(a),std(a)));
            end
        end
        suptitle([st ' r-angle by room(rm#) rscore> ' n2(rst)]);
        
        
        % SCATTER Figure pre vs post pre vs dur rangle
        rst=0.3;
        x=1;  pr = [200-10*x 10 800 700]; 
        figure(900+x); clf; %set(gcf,'position',pr);
        t = [cellsn.(s{x})]; ts = [t.(vs)]; tl=arrayfun(@(z) len(z.st),t); tis = ts>rst&tl>100;
        t = [cellsn.(s{1})]; tab=rad2deg([t(tis).(va)])'; %add tl's of all sessions??
        t = [cellsn.(s{2})]; tam=rad2deg([t(tis).(va)])';
        t = [cellsn.(s{3})]; taa=rad2deg([t(tis).(va)])';
        subplot(121); plot(tab,tam,'o'); axis equal; axis square;
        title('r-angle'); xlabel('pre'); ylabel('dur');
        subplot(122); plot(tab(taa~=0),taa(taa~=0),'o'); axis equal; axis square;
        title('r-angle'); xlabel('pre'); ylabel('post');
        suptitle(['r-score>' n2(rst) ' ' tss{2}]);        
        
        
        
        
        
        % WINDOW ANGLE Figure
        rst=0.3; for x=1:3; pr = [200-10*x 10 800 700];
        t = [cellsn.(s{x})]; tl=arrayfun(@(z) len(z.st),t); 
        ts = [t.(vs)];  tis = ts>rst & tl>100; t=t(tis); %ta = rad2deg([t.(va)]);
        [ sum(tis)]
        figure(110+x); clf; %set(gcf,'position',pr);
        b = 0; y = 5:10:35;
        for z = y
            e=[];
            for a = t
                %e(end+1)=a.pt(end)/60;
                c=windowsesh(a,0,0, z*60, (z+10)*60);
                if ~isempty(c)
                    %[e(end+1,1), e(end,2)]= maxx(rand(1e4,1));
                    %[e(end+1,1), e(end,2)]=rayleigh_score(c.pt,c.st,c.hd); %end+1 did not work????
                    f = size(e,1)+1;
                    [e(f,1), e(f,2)] = rayleigh_score(c.pt,c.st,c.hd); %end+1 did not work????
                    
                end
            end
            
            b=b+1; subplot( ceil(len(y)/2),ceil(len(y)/2),b);
            if ~isempty(e);fpp(s{x},z,e(:,2));end
        end
        suptitle(sprintf('%s: r-angle histogram by window r-rcore>%.1f',tss{x},rst));
        end
        
        function fpp(s,z,a)
            a = rad2deg(a); cla;
            histogram(a,-180:30:180); xlim([-180 180]); ylim([0 60]); xlabel('r angle');
            title(sprintf('%d:%dmin n=%d a=%.2f +- %.2f',z,z+10,len(a),mean(a),std(a)));
        end
        
        %sort(unique(e))'
        %stats(e)
        
        
        
        %works for finding groups
        %a = [clls(ismember([clls.ind],gh1)).midall];
        hdgs = num2cell(findgroups({hdclls.id}, {hdclls.date}));
        [hdclls.group] = hdgs{:};
        gs = unique( str2double(strcat({cellsn.id}, {cellsn.date})) );
        atop = arrayfun(@(x) cellsn(str2double(strcat({cellsn.id}, {cellsn.date}))==x),gs,'uni',0);
        ai= cellfun (@(x) [x.ind],atop,'uni',0)';
        ctop = arrayfun(@(x) hdclls([hdclls.group]==x),unique([hdclls.group]),'uni',0);%cell per group
        ci= cellfun (@(x) [x.ind],ctop,'uni',0)';
        
        % {
        %ALL GROUPS
        figure(201);clf; t = rand(40+3,3); %jet(40+3);
        
        subplot(121);cla; hold on; v = 'rayleigh_angle';
        set(gca, 'ColorOrder', t),
        a = cellfun (@(x) [x.(s{1})],atop,'uni',0)';
        a = cellfun (@(x,q) [q;rad2deg([x.(v)])]',a,ai,'uni',0);
        cellfun(@(x) plot(x(:,1),x(:,2),'.','markersize',25), a);
        %b = [  clls.(s)]; bi = [clls.ind];   b = rad2deg([b.(v)]);plot(bi,b,'o');
        c = cellfun (@(x) [x.(s{1})],ctop,'uni',0)';
        c = cellfun (@(x,q) [q;rad2deg([x.(v)])]',c,ci,'uni',0);
        cellfun(@(x) plot(x(:,1),x(:,2),'o','markersize',25), c);
        %legend({'all';'cohort';'hd g1'; 'hd g2'; 'hd g3'});
        xlabel('cell id'); title(sprintf('%s %s',tss{1},v)); axis tight;
        
        subplot(122); hold on; v = 'rayleigh_angle';
        set(gca, 'ColorOrder', t),
        a = cellfun (@(x) [x.(s{2})],atop,'uni',0)';
        a = cellfun (@(x,q) [q;rad2deg([x.(v)])]',a,ai,'uni',0);
        cellfun(@(x) plot(x(:,1),x(:,2),'.','markersize',25), a);
        %b = [  clls.(s)]; bi = [clls.ind];   b = rad2deg([b.(v)]);plot(bi,b,'o');
        c = cellfun (@(x) [x.(s{2})],ctop,'uni',0)';
        c = cellfun (@(x,q) [q;rad2deg([x.(v)])]',c,ci,'uni',0);
        cellfun(@(x) plot(x(:,1),x(:,2),'o','markersize',25), c);
        %legend({'all';'cohort';'hd g1'; 'hd g2'; 'hd g3'});
        xlabel('cell id'); title(sprintf('%s %s',tss{2},v)); axis tight;
        
        
        %COMPLETE HD GROUPS
        figure(200);clf;
        
        subplot(221); cla; hold on; v = 'rayleigh_angle';
        a = [cellsn.(s{1})];                     a = rad2deg([a.(v)]);plot(a,'.');
        b = [  clls.(s{1})]; bi = [clls.ind];   b = rad2deg([b.(v)]);plot(bi,b,'x');
        c = cellfun (@(x) [x.(s{1})],ctop,'uni',0)';
        c = cellfun (@(x,q) [q;rad2deg([x.(v)])]',c,ci,'uni',0);
        cellfun(@(x) plot(x(:,1),x(:,2),'o'), c);
        legend({'all';'cohort';'hd g1'; 'hd g2'; 'hd g3'});
        xlabel('cell id'); title(sprintf('%s %s',tss{1},v)); axis tight;
        
        subplot(222); hold on; v = 'rayleigh_score';
        a = [cellsn.(s{1})];                     a = [a.(v)];plot(a,'.');
        b = [  clls.(s{1})]; bi = [clls.ind];   b = [b.(v)];plot(bi,b,'x');
        c = cellfun (@(x) [x.(s{1})],ctop,'uni',0)';
        c = cellfun (@(x,q) [q;[x.(v)]]',c,ci,'uni',0);
        cellfun(@(x) plot(x(:,1),x(:,2),'o'), c);
        xlabel('cell id'); title(sprintf('%s %s',tss{1},v)); axis tight;
        
        subplot(223); hold on; v = 'rayleigh_angle';
        a = [cellsn.(s{2})];                     a = rad2deg([a.(v)]);plot(a,'.');
        b = [  clls.(s{2})]; bi = [clls.ind];   b = rad2deg([b.(v)]);plot(bi,b,'x');
        c = cellfun (@(x) [x.(s{2})],ctop,'uni',0)';
        c = cellfun (@(x,q) [q;rad2deg([x.(v)])]',c,ci,'uni',0);
        cellfun(@(x) plot(x(:,1),x(:,2),'o'), c);
        xlabel('cell id'); title(sprintf('%s %s',tss{2},v)); axis tight;
        
        subplot(224); hold on; v = 'rayleigh_score';
        a = [cellsn.(s{2})];                     a = [a.(v)];plot(a,'.');
        b = [  clls.(s{2})]; bi = [clls.ind];   b = [b.(v)];plot(bi,b,'x');
        c = cellfun (@(x) [x.(s{2})],ctop,'uni',0)';
        c = cellfun (@(x,q) [q;[x.(v)]]',c,ci,'uni',0);
        cellfun(@(x) plot(x(:,1),x(:,2),'o'), c);
        xlabel('cell id'); title(sprintf('%s %s',tss{2},v)); axis tight;
        
        %{
        %plot hd per cell %LARGE PLOTS
        x=1;
        rst=0.3; str.x = ''; str.y = '';
        t = [cellsn.(s{x})]; tl=arrayfun(@(z) len(z.st),t)'; ts = [t.(vs)]'; t=t(ts>rst&tl>100); 
        [~, y] = sort(wrapTo360(rad2deg([t.(va)]))); t=t(y);
        n=sort(arrayfun(@(f) len(f.st),t ))';
        sprintf('%s n=%d a=%.2f +- %.2f',s{x},len(n),mean(n),std(n))
        
        figure(400+x); clf; pr = [100 10 1600 1000]; set(gcf,'position',pr);
        for i=1:len(t)
            subplot(ceil(sqrt(len(t))),ceil(sqrt(len(t))),i)
            z = t(i);
            str.t = sprintf('%d%c',wrapTo360(round(rad2deg(z.(va)))),char(176));
            plotHD(gca,z,str); axis square;
        end
        suptitle([s{x} ': rayleigh angle [0=North 90=West] by cell rst> ' n2(rst)]);
        
        figure(4000+x); clf; pr = [100 10 1600 1000]; set(gcf,'position',pr);
        for i=1:len(t)
            subplot(ceil(sqrt(len(t))),ceil(sqrt(len(t))),i)
            z = t(i);
            [z.(va), z.hd, ~, ~] = movingDirection(z,10);
            %c.si = discretize(c.st, [-Inf; mean([c.pt(2:end) c.pt(1:end-1)],2); +Inf]);
            %rd = histcounts(wrapTo360(rad2deg(c.hd(c.si))),bins)./...
            str.t = sprintf('%d%c',wrapTo360(round(rad2deg(z.(va)))),char(176));
            plotHD(gca,z,str);
        end
        suptitle([s{x} ': moving direction angle [0=North 90=West] by cell rst> ' n2(rst)]);
        
        x=1; t = [cellsn.(s{x})];
        figure(300+x); clf; pr = [100 10 1600 1000]; set(gcf,'position',pr);
        for i=1:len(t)
            subplot(ceil(sqrt(len(t))),ceil(sqrt(len(t))),i)
            z = t(i);
            z.si = discretize(z.st, [-Inf; mean([z.pt(2:end) z.pt(1:end-1)],2); +Inf]);
            z.shd=z.hd(z.si);
            [zz.(va), zz.hd, zz.shd, zz.(vs)] = movingDirection(z,10);
            %rd = histcounts(wrapTo360(rad2deg(z.shd),bins)./...
            a=rad2deg(z.hd-zz.hd);histogram(a,-180:1:180);xlim([-180 180]);
            title(sprintf('rd%d%c md%d%c ',round(rad2deg(z.(va)-zz.(va))),char(176),...
                round(rad2deg(mean(z.hd-zz.hd))),char(176)) );
        end
        suptitle([s{x} ': HD-MD diff angle']);
        
%         unique(a)
%         gh1 = [100,101,103,104,105,111]; % #29
%         gh2 = [112,114,115,116,117,118,119];%ci = [112, 114, 117, 118, 119]; #30
%         gh3 = [282,283,287]; %#15
        
        %}
        
    end %hadangle()

    function fieldAngle()
        a=[],b=[],c=[],d=[],e=[],f=[],m=[],n=[],x=[],y=[],z=[];
        
        
        %pre vs during scatter, and pre vs post, first 20
        cang=[];
        for j=1:len(pairs)
            j
            c = aclls(pairs(j,1)); d = aclls(pairs(j,2));
            [~, a]=gridscore2(xcorr2nan(c.before.rm,d.before.rm),2);
            [~, b]=gridscore2(xcorr2nan(c.midall.rm,d.midall.rm),2);
            [~, e]=gridscore2(xcorr2nan(c.after.rm, d.after.rm) ,2);
            %[wrapTo180(a) wrapTo180(b) wrapTo180(e)]
            cang(j,:) =  [wrapTo180(a) wrapTo180(b) wrapTo180(e)];
        end
        
        all=wrapTo180(cang(:,2)-cang(:,1));
        all2=wrapTo180(cang(:,3)-cang(:,1));
        
        figure(232212); clf
        subplot(222)
        hold off;
        bar(sort(abs(all))); b= gca(); %b.YLim=[0 180];
        b.XLabel.String = 'pair'; b.YLabel.String='angle';
        b.Title.String = 'difference in angle pre vs dur';
        subplot(224)
        hold off;
        bar(sort(abs(all2))); b= gca(); %b.YLim=[0 180];
        b.XLabel.String = 'pair'; b.YLabel.String='angle';
        b.Title.String = 'difference in angle pre vs post';
        subplot(221); hold on;  axis square;
        d1=cang(:,1);d2=cang(:,2); plot(d1,d2,'.');
        b= gca(); b.XLabel.String = 'pre'; b.YLabel.String='dur';
        b.Title.String = 'angle pre vs dur';
        f1 = fit(d1,d2,'poly1');%plot(d1,f1(d1),'-'); axis('tight'); axis square;
        text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
        subplot(223); hold on;  axis square;
        d1=cang(:,1);d2=cang(:,3); plot(d1,d2,'.');
        b= gca(); b.XLabel.String = 'pre'; b.YLabel.String='post';
        b.Title.String = 'angle pre vs post';
        f1 = fit(d1,d2,'poly1');%plot(d1,f1(d1),'-'); axis('tight'); axis square;
        text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f1.p1,3), round(ccof(d1,d2),3)),'Units','normalized');
        
        suptitle('Cell pair cross-correlation: angle of first firing field');
        
        %(xA-xB,yA-yB) for |tB-tA|<1s,
        %'xi,yi' the position of the animal at times 'ti' of neuron 'i' spikes.
    end %fieldAngle()

    function misc %Questions
        a=[];aa=a;b=[];bb=a;c=[];cc=a;d=[];e=[];f=[];g=[];h=[];m=[];q=[];t=[];w=[];z=[];
        % misc
        
        % Average Angle
        % vector sum of angle: add all sins and all coss of each angle, arctan of all; or mean sins / mean cos
        x=2; t = [cellsn.(s{x})]; tl=arrayfun(@(z) len(z.st),t);   
        ta = [t.(va)]; ts = [t.(vs)];  tis = ts>rst&tl>=100; ta = ta(tis);
        rad2deg(atan(sum(cos(ta))/sum(sin(ta))))
        
        
        
        %was significance stable
        n = max(ptsbma(:,1));
        bar([sum(ptsbma(:,1)<= n*0.01)/length(pairs) + sum(n*0.99 <= ptsbma(:,1))/length(pairs);...
            sum(ptsbma(:,2)<= n*0.01)/length(pairs) + sum(n*0.99 <= ptsbma(:,2))/length(pairs);...
            sum(ptsbma(:,3)<= n*0.01)/length(pairs) + sum(n*0.99 <= ptsbma(:,3))/length(pairs);]*99,'stacked')
        a = 71/99; b = 45/99;(a-b)^2/a ; [h p]=chi2gof([a b]);
        a = 71; b = 45;(a-b)^2/a ;
        
        a = 58/704; b = 32/704;(a-b)^2/a ; [h p]=chi2gof([a b]);
        a = 58; b = 32;(a-b)^2/a ;
        
        
        %71/99 %45/99
        
        %were grid->hd mixed with grid-/->hd
        for i = 1:len(ctop)
            c = ctop{i};
            b=[c.before]; c=[c.midall];
            [[b.gridscore]>0.3; [c.rayleigh_score]]
        end
        t = [29 30 15];
        for i = 1:3
            c = ctop{i};%atop{t(i)};%ctop{i}; %for hd cohort only
            b=[c.before]; c=[c.midall];
            [[b.gridscore]>0.3; [c.rayleigh_score]]
        end
        
        %number of animals;
        len(unique(str2double({clls.id})))
        %Wilcoxon
        [p,h,wstats]=signrank(ctsbma(:,1),ctsbma(:,2)) %pre mid
        
        [p,h,wstats]=signrank(ctsbma(ctsbma(:,3)~=0,1),ctsbma(ctsbma(:,3)~=0,3)) %pre post
        %h=1 - rejection of the null hypothesis that x-y comes from a distribuation
        %whose median [difference] is 0 (they dont come from the same distribution).
        
        aclls = [     cells{:}];
        %gs = num2cell(findgroups( strcat({aclls.id}, {aclls.date}) ));
        gs = unique( str2double(strcat({aclls.id}, {aclls.date})) );
        atop = arrayfun(@(x) aclls(str2double(strcat({aclls.id}, {aclls.date}))==x),gs,'uni',0);
        ai= cellfun (@(x) [x.ind],atop,'uni',0)';
        
        
        
        
        % EXAMPLES
        fun = @(x) celldata(x).CN == 4 % useful for complicated fields
        tf2 = arrayfun(fun, 1:numel(celldata))
        
        
        A=[3 4 5 6 7];
        B=[6 4 7];
        [sharedvals,idx] = intersect(A,B,'stable'); ismember(A,B)
        
        %a = [clls(ismember([clls.ind],gh1)).before];
        
        V = [2 0 5 7 0 14 0 17 19 20];
        idx = [find(V == 0) length(V)+1];                   % Find Zeros & End Indices
        didx = diff([0 idx])-1;                             % Lengths Of Nonzero Segments
        V0 = V;                                             % Create ‘V0’ (Duplicate Of ‘V’)
        V0(V0 == 0) = [];                                   % Delete Zeros
        C = mat2cell(V0, 1, didx);                          % Create Cell Array
        Result = C(cellfun(@(x) size(x,2)~=0, C));
        
        A = [10,3,1;10,4,2;10,1,1;10,2,1;10,1,3;10,2,4;10,2,5];
        U = unique(A(:,2));
        C = arrayfun(@(x)A(A(:,2)==x,:),U,'uni',0); %'uni' can return any value (returns cells)
        C{:}
        
    end


end %END FILE

%{


         %REMOVE 0 from middle
        %{ 
        %1or 2
        b= 2;  mm = []; dxywinr=dxywinro; i=len(dxywinr{1,1}.(s{1})); %del=51-b:51+b;
        [y x] = meshgrid(1:i, 1:i); del=(x-51).^2+(y-51).^2<b^2; %rad=b;
        for twin = twins
            for i = 1:len(dxywinr)
                for z = s
                    %t=dxywinr{i,twin}.(z{:}); t(del,del)=nan; dxywinr{i,twin}.(z{:})=t;
                    dxywinr{i,twin}.(z{:})(del)=nan;
                    t=dxywinr{i,twin}.(z{:});
                    [~, t] = max(t(:));mm(end+1)=t;
                    %t=std(t(~isnan(t(:))));mm(end+1)=t;
                    %t = isoutlier(t(~isnan(t(:)))); mm(end+1)=sum(t);
                end
            end
        end;mm=mm';
        [m f]=mode(mm);%[mean(mm) std(mm)]
        [y x] =ind2sub(size(dxywinr{i,twin}.midall),m);
        {'ind freq unq len' x y f len(unique(mm)) len(mm)}%/len(mm)*100        
        [b a]=arrayfun(@(x) ind2sub([100 100],x), mm); a=a';b=b';
        [mean(a) std(a) mean(b) std(b)]
        %}





        %c = cellfun (@(x) [x.midall],c,'uni',0)';
        %cs = cellfun (@(x,q) [q;       ([x.rayleigh_score])]',c,ci,'uni',0);
        
        
        rsp=[];rs = [];
        for i=1:length(pairs)
            rsp(i,:) = [cells{pairs(i,1)}.before.rayleigh_score,...
                cells{pairs(i,2)}.before.rayleigh_score,...
                cells{pairs(i,1)}.midall.rayleigh_score,...
                cells{pairs(i,2)}.midall.rayleigh_score];
        end
        for i=1:length(cels)
            rs(i,:)  = [cells{cels(i)}.before.rayleigh_score,...
                cells{cels(i)}.midall.rayleigh_score];
        end
        ccl2 = rs(:,2)>0.4; ccl1 = ~ccl2; chd = cels(ccl2);
        cl2 = rsp(:,3)>0.4 & rsp(:,4)>0.4; cl1 = ~cl2;
        
        plot(rs(:,1),rs(:,2),'b.','markersize',mfs);
        hold on; plot(rs(rs(:,2)>0.5,1),rs(rs(:,2)>0.5,2),'r.','markersize',mfs);
        
        a = {};
        for i=1:length(cels)
            a{end+1}=cells{cels(i)}.id;
        end
        unique(a)
        
        k2 = cels(ccl2);
        a = {}; c = struct;
        for i=1:length(k2)
            s = strcat('s',cells{k2(i)}.date, cells{k2(i)}.id);
            if isfield(c, s)
                c.(s) = [c.(s) k2(i)];
            else
                c.(s) = k2(i);
            end
            a{end+1}=strcat(s);
        end
        b = struct2cell(c);
        


        %% drift
         while tsi <= len(b.st)
                            if   b.st(tsi) + twin < t
                                tssi = tsi;
                            elseif  t < b.st(tsi)-twin
                                tsi=len(b.st);
                            else
                                %if b.st(k)-twin <= t && t <= b.st(k)+twin
                                dx =  b.sx(tsi) - a.sx(i);
                                dy =  b.sy(tsi) - a.sy(i);
                                dxys(end+1,:) = [dx dy];
                            end
                            tsi=tsi+1;
                        end
        while tpi <= len(b.pt)
                            if   b.pt(tpi) + twin < t
                                tppi = tpi;
                            elseif  t < b.pt(tpi)-twin
                                tpi=len(b.pt);
                            else
                                %if b.st(k)-twin <= t && t <= b.st(k)+twin
                                dx =  b.px(tpi) - a.sx(i);
                                dy =  b.py(tpi) - a.sy(i);
                                dxyt(end+1,:) = [dx dy];
                            end
                            tpi=tpi+1;
                        end
        
        
        
        
        for h = 1:len(twins)
            twin = twins(h)
            for j=1:len(pairs)
                tic
                %BEFORE
                a = cells{pairs(j,1)}.before; b = cells{pairs(j,2)}.before;
                dxy = [];
                ksi = 1;
                for i = 1:len(a.st)
                    t = a.st(i);
                    k = ksi;
                    while k <= len(b.st)
                        if   b.st(k) + twin < t
                            ksi = k;
                        elseif  t < b.st(k)-twin
                            k=len(b.st);
                        else
                            %if b.st(k)-twin <= t && t <= b.st(k)+twin
                            dx =  b.sx(k) - a.sx(i);
                            dy =  b.sy(k) - a.sy(i);
                            dxy(end+1,:) = [dx dy];
                        end
                        k=k+1;
                    end
                end
                dxywin{j,twin}.b = dxy;
                %DURING
                a = cells{pairs(j,1)}.midall; b = cells{pairs(j,2)}.midall;
                dxy = [];
                ksi = 1;
                for i = 1:len(a.st)
                    t = a.st(i);
                    k = ksi;
                    while k <= len(b.st)
                        if   b.st(k) + twin < t
                            ksi = k;
                        elseif  t < b.st(k)-twin
                            k=len(b.st);
                        else
                            dx =  b.sx(k) - a.sx(i);
                            dy =  b.sy(k) - a.sy(i);
                            dxy(end+1,:) = [dx dy];
                        end
                        k=k+1;
                    end
                end
                dxywin{j,twin}.m = dxy;
                toc
            end
        end
        
        dbmg = {};
        for h = 1:len(twins)
            twin = twins(h);
            t = [];
            for i = 1:len(pairs)
                [i twin]
                tic
                nbins = 50;
                a = dxywin{i,twin};
                px = a.b(:,1);py = a.b(:,2); %px = px+max(px); py = py+max(py);
                mx = max(px);
                my = max(py);
                pxi = discretize(px, min(px):range(px)/nbins:max(px));
                pyi = discretize(py, min(py):range(py)/nbins:max(py));
                rmt = accumarray([pyi, pxi], 1, [nbins nbins]);
                gb = gridscore2(xcorr2(rmt),2);
                nbins = 50;
                a = dxywin{i,twin};
                px = a.m(:,1);py = a.m(:,2); %px = px+max(px); py = py+max(py);
                mx = max(px);
                my = max(py);
                pxi = discretize(px, min(px):range(px)/nbins:max(px));
                pyi = discretize(py, min(py):range(py)/nbins:max(py));
                rmt = accumarray([pyi, pxi], 1, [nbins nbins]);
                gm = gridscore2(xcorr2(rmt),2);
                toc
                t(i,:)= [gb gm];
            end
            dbmg{twin} = t;
            h
        end
        
        
        
        
        dbmc = {};
        for h = 1:len(twins)
            twin = twins(h);
            t = [];
            for i = 1:len(pairs)
                %ab = cells{pairs(i,1)}.before.rm; bb = cells{pairs(i,2)}.before.rm;dbmca(i) = ccof(ab-bb,am-bm);
                t(i)= driftccof(dxywin{i,twin}.b, dxywin{i,twin}.m);
            end
            dbmc{twin} = t;
            h
        end
        figure(27);
        plot(1:99,[sort(dbmc{1});sort(dbmc{3});sort(dbmc{5});sort(dbmc{10})],'.');
        title('corr pre dur pair drift');
        legend({'1s';'3s';'5s';'10s'});
        
        figure(28); hold on;
        for h = 1:len(twins)
            subplot(2,2,h)
            twin = twins(h);
            a = dbmg{twin};
            plot(a(:,1),a(:,2),'.');
            title(sprintf('drift %ds',twin)); axis square; axis equal;
            xlabel('gridscore pre'); ylabel('gridscore dur'); mm=slim(gca);
            xlim(mm);ylim(mm)
        end
        suptitle('gridscore of 2D dxdy histogram by drift window');
        
        % graphics imagesc
        for h = 1:len(twins)
            twin = twins(h); g = dbmg{twin};
            figure(100+twin)
            for i = 1:20%(pairs)
                nbins = 50;
                a = dxywin{i,twin}.b;
                px = a(:,1);py = a(:,2); %px = px+max(px); py = py+max(py);
                pxi = discretize(px, min(px):range(px)/nbins:max(px));
                pyi = discretize(py, min(py):range(py)/nbins:max(py));
                rmt = accumarray([pyi, pxi], 1, [nbins nbins]);
                subplot(4,10,2*i-1)
                imagesc(imgaussfilt(rmt, 2,'FilterDomain','spatial'));
                colormap jet; axis off; axis square;title(sprintf('p%dpre:%.2f',i,g(i,1)));
                a = dxywin{i,twin}.m;
                px = a(:,1);py = a(:,2); %px = px+max(px); py = py+max(py);
                pxi = discretize(px, min(px):range(px)/nbins:max(px));
                pyi = discretize(py, min(py):range(py)/nbins:max(py));
                rmt = accumarray([pyi, pxi], 1, [nbins nbins]);
                subplot(4,10,2*i)
                imagesc(imgaussfilt(rmt, 2,'FilterDomain','spatial'));
                colormap jet; axis off; axis square;title(sprintf('p%ddur:%.2f',i,g(i,2)));
            end
            suptitle(sprintf('dxdy hist and gridscore by pair %ds',twin));
        end
        
        
        figure(26)
        for i = 1:15%(pairs)
            nbins = 50;
            a = dxym10{i};
            px = a(:,1);py = a(:,2); %px = px+max(px); py = py+max(py);
            mx = max(px);
            my = max(py);
            pxi = discretize(px, min(px):range(px)/nbins:max(px));
            pyi = discretize(py, min(py):range(py)/nbins:max(py));
            rmt = accumarray([pyi, pxi], 1, [nbins nbins]);
            subplot(3,5,i)
            imagesc(rmt); colormap jet; axis off; axis square;title(sprintf('pair%d',i));
        end
        suptitle('dxdy dur 10s');
        
        figure(23)
        for i = 1:15%(pairs)
            nbins = 50;
            a = dxyb{i};
            px = a(:,1);py = a(:,2); %px = px+max(px); py = py+max(py);
            mx = max(px);
            my = max(py);
            pxi = discretize(px, min(px):range(px)/nbins:max(px));
            pyi = discretize(py, min(py):range(py)/nbins:max(py));
            rmt = accumarray([pyi, pxi], 1, [nbins nbins]);
            subplot(3,5,i)
            imagesc(rmt); colormap jet; axis off; axis square;title(sprintf('pair%d',i));
        end
        suptitle('dxdy pre 1s');
        
        figure(24)
        for i = 1:15%(pairs)
            nbins = 50;
            a = dxym{i};
            px = a(:,1);py = a(:,2); %px = px+max(px); py = py+max(py);
            mx = max(px);
            my = max(py);
            pxi = discretize(px, min(px):range(px)/nbins:max(px));
            pyi = discretize(py, min(py):range(py)/nbins:max(py));
            rmt = accumarray([pyi, pxi], 1, [nbins nbins]);
            subplot(3,5,i);
            imagesc(rmt); colormap jet; axis off; axis square; title(sprintf('pair%d',i));
        end
        suptitle('dxdy dur 1s');
        
        for j=1:len(cels)
            c = cells{cels(j)}.midall.pt;
            a(j)=c(1)/60;
        end
        a = unique(a)';

       
        %gridwin={};
        %histogram scatters
        dxyr=dxywinrd; z='dmet';%dxyr=dxywinrn; z='nmet';
        m=1;w=2; mm=['corr ' ts{m} ' ' ts{w}]; g=[];%compare against
        figure(3000+w); clf; pr = [200 10 800 600]; set(gcf,'position',pr);
        for h =1:len(twins)
            twin=twins(h);
            %g={}; %scatter
            %for w = 2:3 %scatter
                a=[];
                for i = 1:len(dxyr)
                tic; [twin i]
                %corr sesh vs sesh for dxdy
                t=corr(dxyr{i,twin}.(s{m})(:),dxyr{i,twin}.(s{w})(:),'rows','complete');
                %gridscore sesh sesh dxdy
                %to make e=dxyr{i,twin}.(s{1}); f=dxyr{i,twin}.(s{2}); d=dxyr{i,twin}.(s{3}); 
                %   e=xcorr2nan(e,e); f=xcorr2nan(f,f);d=xcorr2nan(d,d);autowin{i,twin}.before=e; 
                %   autowin{i,twin}.midall=f; autowin{i,twin}.after=d;
                %t=[gridscorenan(e,5) gridscorenan(f,5)]
                a(end+1,:)=t; 
                toc
                end
            g(:,h)=a; %g{h}=a(~isnan(a));%
            %g{w}=a; %scatter
            %gridwin{twin}=a;
            %end
           %subplot(2,3,h);
            %histogram
           %histogram(a,0:0.02:0.5);
           %title(['dxdy drift win ' n2(twin) 's']);
           %xlabel(mm); axis square;
            %scatter
%             y=g{2}';x=g{3}';t=~isnan(y)&~isnan(x);x=x(t);y=y(t); scatter(x,y); hold on;
%             f = fit(x, y,'poly1');plot(x,f(x),'-'); axis('tight'); axis square;
%             text(0.1,0.9,sprintf('a=%.3f r=%.3f',round(f.p1,3), round(ccof(x,y),3)),'Units','normalized');
%             title([n2(twin) 's']); xlabel('pre vs post'); ylabel('pre vs dur'); 
%             e=slim(gca);xlim(e);ylim(e);
        end
        %histogram
        g(any(isnan(b), 2), :) = []; 
        clf;violin(g,'facecolor',summer(5));ylabel(mm); xlabel('window (sec)'); 
        xticks(1:len(twins));xticklabels(twins);
        t=[mm ' drift dxdy ' z];
        %t=[mm ' drift dxdy ' z];
        suptitle(t);
        
        saveas(gcf,['./figs/violin ' t '.png']);
%         saveas(gcf,['./figs/hist ' t '.png']);
        %scatter
%         t=['correlation of dxdy by time window ' z];
%         suptitle(t);
%         saveas(gcf,['./figs/scatter 't '.png']);     




b='';
for i = 1:length(files) %CHECK 36
    a=cells{i};
    if ~isequal(b,strcat(a.id,a.date))
        b=strcat(a.id,a.date);a=a.midall.pt; c(end+1,:)=[min(a),max(a)];
        %fprintf('%d %.1f %.1f \n',i, min(a)/60, max(a)/60);
    end
end
c=[]; b='';
for j = 1:length(cels) %CHECK 36
    i=cels(j);
    a=cells{i};
    if ~isequal(b,strcat(a.id,a.date))
        b=strcat(a.id,a.date);a=a.midall.pt; c(end+1,:)=[min(a),max(a)];
        fprintf('%d %.1f %.1f \n',i, min(a)/60, max(a)/60);
    end
end
    

%}



