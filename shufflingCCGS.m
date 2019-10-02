function shufflingCCGS(cells)
    pairs = [2,7;11,12;11,13;11,14;12,13;12,14;13,14;22,24;22,28;24,28;34,35;34,37;34,39;34,41;35,37;35,39;35,41;37,39;37,41;39,41;51,53;61,62;61,66;61,68;62,66;62,68;66,68;77,78;77,79;78,79;80,81;82,86;100,101;100,103;100,104;100,105;100,111;101,103;101,104;101,105;101,111;103,104;103,105;103,111;104,105;104,111;105,111;112,114;112,115;112,116;112,117;112,118;112,119;114,115;114,116;114,117;114,118;114,119;115,116;115,117;115,118;115,119;116,117;116,118;116,119;117,118;117,119;118,119;171,173;171,174;171,176;171,177;173,174;173,176;173,177;174,176;174,177;176,177;178,179;178,180;178,182;178,184;179,180;179,182;179,184;180,182;180,184;182,184;188,189;199,201;199,209;201,209;230,231;274,276;274,278;276,278;282,283;282,287;283,287];
    %SHUFFLING Spatial + time
    
    par=[]; par.n=10; par.movmean=25; par.nb = 1000; 
    
    t=zeros(len(pairs),par.n);pgcb = t; pgcm = t; pgca = t;
    
    for i = 1:length(pairs)
        tic
        disp(i)
        c1 = cells{pairs(i,1)};
        c2 = cells{pairs(i,2)};
        %space
        pgcb(i,:) = shuffleCCGS(c1.before, c2.before,par); %[pairs(i,:), ]
        pgcm(i,:) = shuffleCCGS(c1.midall, c2.midall,par);
        if(length(c1.after.pt)>1 && length(c2.after.pt)>1)
            pgca(i,:) = shuffleCCGS(c1.after, c2.after,par);
        else;pgca(i,:) = zeros(1,par.n);end
         toc
    end
    
    sum(pgcm>.1)
    sum(pgcb>.1)
    
    %return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tt = {ptcb,ptcm,ptca,pscb,pscm,psca,pgcb,pgcm,pgca};
    ctsbma = zeros(len(pairs),len(tt));
    ptsbma = zeros(len(pairs),len(tt));
    lp = len(pairs);
    for ii = 1: len(tt)
        ctsbma(:,ii) = tt{ii}(:,1);
        for i = 1:lp
            t = (tt{ii}(i,:));
            [~,I] = sort(abs(t),'descend');
            if ii > 3 %for grid scores
                [~,I] = sort(t,'descend');
            end
            pp = find(I==1); %index of non shuffled
            ptsbma(i,ii) = pp;
            if t(1)==0 || t(1)==-2
                ptsbma(i,ii) = len(t);
            end
        end
    end
    pp = ptsbma;
    pp(pp==1)=par.n*0.01;%for runs smaller than 100
    pp(pp>par.n*0.01)=0; %indices greater than 1
    pp(pp>0)=1;
    pptsbma = logical(pp);
    sum(ptsbma(:,:)<=par.n*0.01|ptsbma(:,:)==1)/length(pairs)*100    
    sum(pptsbma)/length(pairs)*100
    
    sum(ctsbma>0)
    
    return  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
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

