% window drift
% f3s4
%based on 'Revisopm.m'
%load('.\data\fs8dxywin');
%load('.\data\fs8_v2_dxywin');

function fs8(cellsn,pairs,dxywinrdshufpairs)

    %load('C:\Noam\Data\muscimol\aclls15min','aclls'); %personal
%   load('.\\data\\dxdyrate','dxywinrd','dxywinrn','gridwin','autowin');
    %load('.\\data\\dxdyrate15','dxywinrd15'); F %
    nbins = 100; twins = [1 2 3 5 10] ; %window in secs
    
    a = [1 10 12];% 11 [1 15 18]; %CELLS TO SHOW %12 13 16
    
    %dxywinrdshufpairs=makeWins(); %to make dxywinrdshufpairs
    
    s = {'before', 'midall','after'}; ts = {'pre','dur','post'};
    fs = 12;
    gp={'Spacing',5, 'BackgroundColor','w','DividerMarkings','off'};
    ap={'visible','off'};pnt='Parent'; 
    tp={'fontweight','bold','fontsize',fs};
    pp={'BackgroundColor','w','bordertype','none'};
    
    fig = figure(81); clf;
    set(fig,'Color','w', 'Position', [600 0 750 720]);
    gAll = uix.GridFlex(pnt,fig,gp{:});
    gTop = uix.GridFlex(pnt,gAll,gp{:});
    gBot = uix.GridFlex(pnt,gAll,gp{:});
     
    axes(pnt,uix.Panel(pnt,gTop,pp{:}),ap{:});
    text(0,0.2,'A',tp{:});
    gA = uix.GridFlex(pnt,gTop,gp{:});   
    %A - imgsc 
    doff = 16; d = 50-doff:51+doff;
    twin=1; dxyr=dxywinrdshufpairs;
    for id = a
        axes(pnt,uicontainer(pnt,gA));
        x=1; c = dxyr{id,twin}.(s{x}); c=c(d,d);
        fp(c,id,doff,x,'k');
        axes(pnt,uicontainer(pnt,gA));
        x=2; c = dxyr{id,twin}.(s{x}); c=c(d,d);
        fp(c,id,doff,x,'r');
    end
    set(gA,'Widths', zeros(1,len(a)*2) -1 )    
    set(gTop,'Heights', [ 25, -1]);
   
    function fp(c,id,doff,x,cr)
        imgsc(c,2); title(sprintf('p%d 1s %s',id,ts{x}),'color',cr);         
        tmax=round(doff*2/5)*5; 
        tt=slim(gca);%t=floor(t/10)*10;xlim(t);ylim(t);
        tt(3)=tt(end);tt(2)=round(tt(end)/2);%dd=1; %t(1)=t(1)+dd; t(end)=t(end)-dd; 
        xticks(tt);yticks(tt); 
        tl={n2(-tmax);'0';n2(tmax)};xticklabels(tl);yticklabels(tl);
        xlabel('cm');ylabel('cm');
        %axis(gca,'square');%set(gca,'ydir','norm');
    end
    
    axes(pnt,uix.Panel(pnt,gBot,pp{:}),ap{:});
    text(0,0.5,'B',tp{:});
    gB = uix.GridFlex(pnt,gBot,gp{:});    
    
    %SCATTER
    m=1;w=2; mm=['corr ' ts{m} ' ' ts{w}];
    for h =1:len(twins)
        twin=twins(h);
        g={}; 
        for w = 2:3 %X mid after
            a=[];
            for i = 1:len(dxyr)
                %[twin i]
                %corr sesh vs sesh for dxdy
                t=corr(dxyr{i,twin}.(s{m})(:),dxyr{i,twin}.(s{w})(:),'rows','complete');
                a(end+1,:)=t;
            end
            g{w}=a;
        end
        axes(pnt,uicontainer(pnt,gB));
        y=toCol(g{2});x=toCol(g{3});%2 mid, 3 after
        t=~isnan(y)&~isnan(x);x=x(t);y=y(t); 
        scatter(x,y); hold on;
        f = fit(x, y,'poly1');%plot(x,f(x),'-'); 
        axis('tight'); 
        [r p]=ccof(x,y);
        %text(0.1,0.9,sprintf('a=%.2f r=%.2f p=%.2f',round(f.p1,2), round(r,2), round(p,2)),'Units','normalized');
        text(0.1,0.9,sprintf('a=%.2f r=%.2f %s',round(f.p1,2), round(r,2), pstr(p,2)),'Units','normalized');
        title([n2(twin) 's']); xlabel('pre vs post corr'); ylabel('pre vs dur corr'); box on
        e=slim(gca);xlim(e);ylim(e);
        plot(e,f(e),'-'); axis square; 
    end
    %scatter
    t=['correlation of dxdy by time window ' twin];
    %suptitle(t);
    set(gB,'Widths', [-1 -1 -1])
    set(gBot,'Heights', [ 25 -1]);
    set(gAll,'Heights', [ -1 -3]);
    a = findobj(fig,'type','UIContainer');
    for i = 1:len(a)
        a(i).BackgroundColor = 'w';
    end
    %saveas(f,['./figs/scatter ' t '.png']);
    
    
%{
            % % % % graphics imagesc 
        d = 16; d = 50-d:51+d; %b=1;
        pr = [200 10 1800 900];  t=2; e = 2;% f=5;f=ceil(sqrt(f));%f-#pairs %t 
        for twin = 1%twins
            figure(1000+10*twin+t); clf; set(gcf,'position',pr);
            for i = 1:20%(pairs)
                subplot(4,e*5,e*i-e+1);x = 1; %PRE 1 DUR 2
                a = dxywinrdshufpairs{i,twin}.(s{x}); a=a(d,d); %b = gridscore2(a,t);
                imagesc(imgaussfilt(a, t,'FilterDomain','spatial'));
                %imagesc(xcorr2(imgaussfilt(a, t,'FilterDomain','spatial')));
                %a(isnan(a))=-1;%imagesc(a);
                colormap jet; axis off; axis square;title(sprintf('p%d:%s',i,ts{x}));
                
                subplot(4,e*5,e*i); x = 2;%PRE 1 DUR 2
                a = dxywinrdshufpairs{i,twin}.(s{x}); a=a(d,d); % b = gridscore2((a),t);
                imagesc(imgaussfilt(a, t,'FilterDomain','spatial'));
                %imagesc(imgaussfilt(xcorr2(a), t,'FilterDomain','spatial'));
                colormap jet; axis off; axis square;title(sprintf('p%d:%s',i,ts{x}),'color','red');
            end
            %suptitle(sprintf('%s %ds: dxdy-spike/dxdy-pos and gridscore by pair sigma=%d',s{sr},twin,t));
            suptitle(sprintf('%ds: dxdy-spike/dxdy-pos and gridscore by pair sigma=%d',twin,t));
            %saveas(gcf, sprintf('./figs/drift_diff_%ds_sig%d_zoom_dmet.png',twin,t) );
        end
    %}

    %{%

    function dxywinrdshufpairs = makeWins()
     s = {'before', 'midall','after'};
    t=[cellsn.(s{1})];d=cellfun(@(x) max(diff(x)),{t.pt},'uni',1);
    t=[cellsn.(s{2})];e=cellfun(@(x) max(diff(x)),{t.pt},'uni',1); off=max([d e]);%0.04
    dxywinrdshufpairs={};
    for h = 1:len(twins)
        twin = twins(h);
        for j=1:len(pairs)
            ['a ' n2(twin) ' ' n2(j)]
            tic
            for ss= [s(:)']
                a = cellsn(pairs(j,1)).(ss{:}); b = cellsn(pairs(j,2)).(ss{:});
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
                dxywinrdshufpairs{j,twin}.(ss{:}) = dxyr;
            end
            toc;
        end
    end
    end %end makeWin()

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

end
