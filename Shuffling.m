function Shuffling(cellsn,pairs)
%load('Z:\data\noam\muscimol\cells15nan');
%shuffling gridscore
     p.nshuf=1000;p.pval=0.02;p.spval=0.01;  
                                            p.nb=50;p.movm=25;p.asig=2;
     pgb3 = []; pgm3 = []; pga3 = [];
     for i = 1:len(cellsn) 
         i
         c=cellsn(i); tic;                                                          
         pgb3(i,:) = shuffleGridscoreNan(c.before,p.nshuf,p.nb,p.movm,p.asig,p.pval,'before');
         t=pgb3(i,:); [~,ix]=sort(t,'descend'); ix1=find(ix==1);
          %[n2(ix1) ' ' n2(t(1))]
         if ix1 <= p.nshuf*p.pval;%index of non shuffled
             ['DO MIDALL ' n2(ix1) ' ' n2(t(1))]
         pgm3(i,:) = shuffleGridscoreNan(c.midall,p.nshuf,p.nb,p.movm,p.asig,p.pval,'midall');
         else
            pgm3(i,:)=zeros(1,p.nshuf); t=shuffleGridscoreNan(c.midall,p.nshuf,p.nb,p.movm,p.asig,p.pval,'midall');pgm3(i,1)=t(1);
         end
         %pga3(i,:) = shuffleGridscoreNan(c.after, p.nshuf,p.nb,p.movm,p.asig,p.pval,'after' );
         save(['../gshuff ' date '-' ts], 'pgb3','pgm3','pga3','p')
        toc         
     end
 %%%%%%%%%%%%%%%%%%%% grid analysis%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    tt = {pgb3,pgm3};ncls=size(pgb3,1);
   %tt = {pgb,pgm};ncls=size(pgb,1);
    cgbma = [];%zeros(len(cellsn),len(tt));
    pgbma = [];%zeros(len(cellsn),len(tt));
    for ii = 1: len(tt)
        cgbma(:,ii) = tt{ii}(:,1);
        for i = 1:ncls
            t = (tt{ii}(i,:));
            %[~,I]=sort(abs(t),'descend');if ii>6;[~,I]=sort(t,'descend');end
            [~,I]=sort(t,'descend');
            pp = find(I==1); %index of non shuffled
            pgbma(i,ii) = pp;
            if t(1)==0 %THIS IS FOR GRID SCORE CHECK
                pgbma(i,ii) = len(t); %do i need this??? ???SDSDWDEWE
            end
        end
    end
    %VISUALIZE
    ppgbma = pgbma <= p.nshuf*p.spval; 
    sum(ppgbma)/ncls
    clear pp; clear n; clear i; clear ii; clear I; clear t; clear tt; 
    



    %shuffling Spatial + Temporal
    
    par=[]; par.n=1000; par.pval=0.02;par.spval=0.01; par.movmean=25; par.nb = 50; %t=zeros(len(pairs),par.n);   
    pscb3 = []; pscm3 = []; psca3 = []; ptcb3 = []; ptcm3 = []; ptca3 = [];
    
    for i = 1:     len(pairs) %<<<<<<<<<<<<<  
        i
        c1 = cellsn(pairs(i,1));%cells{pairs(i,1)};
        c2 = cellsn(pairs(i,2));
        %space
        tic
disp('sb');pscb3(i,:) = shuffleSpaceCorrelations(c1.before, c2.before,par); 
disp('sm');pscm3(i,:) = shuffleSpaceCorrelations(c1.midall, c2.midall,par);
        if(length(c1.after.pt)>1 && length(c2.after.pt)>1) %after
disp('sa');psca3(i,:) = shuffleSpaceCorrelations(c1.after, c2.after,par);
        else;psca3(i,:) = zeros(1,par.n);end
        toc
        %time
        tic
disp('tb');ptcb3(i,:) = shuffleTimeCorrelations (c1.before, c2.before,par);
disp('tm');ptcm3(i,:) = shuffleTimeCorrelations (c1.midall, c2.midall,par);
        if(length(c1.after.pt)>1 && length(c2.after.pt)>1)%after
disp('ta');ptca3(i,:) = shuffleTimeCorrelations (c1.after, c2.after,par);
        else;ptca3(i,:) = zeros(1,par.n);end 
        toc
        save(['../tsshuff ' date '-' ts], 'pscb3','pscm3','psca3','par')
    end    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %tt = {ptcb3,ptcm3,ptca3,pscb3,pscm3,psca3}%,pgcb,pgcm,pgca};
    tt = {ptcb,ptcm,ptca,pscb,pscm,psca}%,pgcb,pgcm,pgca};
    niter=size(tt{1},1);
    ctsbma = zeros(niter,len(tt));
    ptsbma = zeros(niter,len(tt));
    for ii = 1: len(tt)
        ctsbma(:,ii) = tt{ii}(:,1);
        for i = 1:niter
            t = (tt{ii}(i,:));
            [~,I]=sort(abs(t),'descend'); %if ii>6;[~,I]=sort(t,'descend');end
            %[~,I]=sort(t,'descend');
            pp = find(I==1); %index of non shuffled
            ptsbma(i,ii) = pp;
            if t(1)==0 || t(1)==-2 %THIS IS FOR GRID SCORE CHECK
                %ptsbma(i,ii) = len(t);
            end
        end
    end
    clear i; clear ii; clear I; clear t; clear tt; clear pp;
    %VISUALIZE
    '(+) corr'
    sum(ptsbma<= par.n*par.spval&ctsbma>0)/niter
    '(-) corr'
    sum(ptsbma<= par.n*par.spval&ctsbma<0)/niter
    %pptsbma = pp;
    pptsbma= ptsbma <= par.n*par.spval;
    'corr'
    sum(pptsbma)/niter
    
    
    
    
    unique(q(:,2))
    figure;histogram(q(:,2),100)
    
    
    q = [];
    for i = 1: size(pscb,1);
        t = abs(ptcb(i,:));
        [~, mi] = max(t(2:end));
        q(i,:) = [t(1), mi];
        %[s,I] = sort(abs(t),'descend');
        %pb = round(find(I==1)/length(I),2); %index of non shuffled
        %[pb pm]; pbs(end+1) = pb; pms(end+1) = pm;
    end
    unique(q(:,2))
    figure;histogram(q(:,2),100)
    
    
    %PLOT MEASUREMENTS OF GRID bef dur after
    b = []; m = []; a = [];
    for i = 1:length(pairs)
        i
        c1 = cells{pairs(i,1)};
        c2 = cells{pairs(i,2)};
        b(i,:) = apear(c1.before, c2.before); %[pairs(i,:), ]
        m(i,:) = apear(c1.midall, c2.midall);
        if(length(c1.after.pt)>1 && length(c2.after.pt)>1)
            a(i,:) = apear(c1.after, c2.after);
        end
    end
    
    h = figure(934); set(h,'Position',[-1000 10 1000 1000]);
    colstr = {'time cc smoothed at 30ms';
          'spatial cc at 0,0, no smoothing';
          'maximum gridscore';
          'gridscore of spatial cc';
         };
    stdd=3;
    t = b;              k=[1:len(t)]; k = t(:,1) < stdd*std(t(:,1));  %leaving out corr outliers
    for i = 1:3
        subplot(3,3,i)
        k = t(:,i+1)>-2 & k;
        scatter(t(k,1),t(k,i+1),'.');
        xlabel(colstr(1));ylabel(colstr(i+1));
        title('before');
    end
    t = m;              k=[1:len(t)]; k = t(:,1) < stdd*std(t(:,1));
    for i = 1:3
        subplot(3,3,3+i)
        k = t(:,i+1)>-2 & k;
        scatter(t(k,1),t(k,i+1),'.');
        xlabel(colstr(1));ylabel(colstr(i+1));
        title('during');
    end
    t = a(a(:,1)~=0,:); k=[1:len(t)]; k = t(:,1) < stdd*std(t(:,1));
    for i = 1:3
        subplot(3,3,6+i)
        k = t(:,i+1)>-2 & k;
        scatter(t(k,1),t(k,i+1),'.');
        xlabel(colstr(1));ylabel(colstr(i+1));
        title('after');
    end
    
    ax = findobj(h,'type','axes');
    for i = 1:length(ax)
        %axis(ax(i),'square');
    end
    
    
end
    
    
    
    


