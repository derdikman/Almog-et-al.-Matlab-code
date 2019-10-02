%
% MAKE FIGURE FOR NB SIZE

%uncomment to load
%{ 
load('Z:\\data\\noam\\Data\\muscimol\\noam\\cells_45min_d_patchtraj_rayleigh'); %lab
load('C:\\Noam\\Data\\muscimol\\noam\\cells_Infmin_d_patchtraj_rayleigh'); %personal
load('C:\\Noam\\Data\\muscimol\\noam\\cells_45min_d_patchtraj_rayleigh'); %personal
%}
load('C:\\Noam\\Data\\muscimol\\noam\\cells_45min_d_patchtraj_rayleigh'); %personal
pairs = [2,7;11,12;11,13;11,14;12,13;12,14;13,14;22,24;22,28;24,28;34,35;34,37;34,39;34,41;35,37;35,39;35,41;37,39;37,41;39,41;51,53;61,62;61,66;61,68;62,66;62,68;66,68;77,78;77,79;78,79;80,81;82,86;100,101;100,103;100,104;100,105;100,111;101,103;101,104;101,105;101,111;103,104;103,105;103,111;104,105;104,111;105,111;112,114;112,115;112,116;112,117;112,118;112,119;114,115;114,116;114,117;114,118;114,119;115,116;115,117;115,118;115,119;116,117;116,118;116,119;117,118;117,119;118,119;171,173;171,174;171,176;171,177;173,174;173,176;173,177;174,176;174,177;176,177;178,179;178,180;178,182;178,184;179,180;179,182;179,184;180,182;180,184;182,184;188,189;199,201;199,209;201,209;230,231;274,276;274,278;276,278;282,283;282,287;283,287];
cels = unique([pairs(:,1),pairs(:,2)])';
%room=cellfun(@(x) str2num(x),{db_all.room})';

%rmi = [];
%rmi=[43];rml=logical(ones(len(pairs),1));rml(rmi)=false; 
clls = [cells{cels}];
%[groups ~] = findSimultaneouslyRecordedCells(cells);
load('.\\data\\1000Nspacetimeshuffle25msooth'  );
%load('.\\data\\ccofNonGroupTimeSpaceNb1000');
%load('.\\data\\figdata');
%load('.\\data\\CCGSshuffling1k'); 
nb = 50;
%load(sprintf('.\\data\\1OutGroupSpaceShuffleNBins%d',nb));
load(sprintf('.\\data\\1kSpaceShuffleNBins%d',nb)); %very similar 1000Nspacetimeshuffle25msooth'??
% pairs(43,:)=[];ctsbma(43,:)=[]; ptsbma(43,:)=[];pptsbma(43,:)=[];cl1(43)=[];cl2(43)=[];zpptsbma(43,:)=[];    
%% COMPARE TO SHUFFLING 
    tt = {ptcb,ptcm,ptca,pscb,pscm,psca}%,pgcb,pgcm,pgca};
    ctsbma = zeros(len(pairs),len(tt));
    ptsbma = zeros(len(pairs),len(tt));
    for ii = 1: len(tt)
        ctsbma(:,ii) = tt{ii}(:,1);
        for i = 1:len(pairs)
            t = (tt{ii}(i,:));
            %[~,I]=sort(abs(t),'descend');if ii>6;[~,I]=sort(t,'descend');end
            [~,I]=sort(t,'descend');
            pp = find(I==1); %index of non shuffled
            ptsbma(i,ii) = pp;
            if t(1)==0 || t(1)==-2 %THIS IS FOR GRID SCORE CHECK
                ptsbma(i,ii) = len(t);
            end
        end
    end
    clear i; clear ii; clear I; clear nb; clear par; clear t; clear tt; clear pp;
    clear ptca; clear ptcb; clear ptcm; clear psca; clear pscb; clear pscm; 
    %VISUALIZE
    n = 1000; show = [1,2,3,4,5,6]; %n = len(ptcb) % ==1000
    pp = ptsbma<= n*0.01;  sum(pp(:,show))/len(pairs) %LARGEST VAL (FIRST DESC)
    pp = n*0.99 <= ptsbma; sum(pp(:,show))/len(pairs) %SMALLEST VAL (LAST DESC)
    pp = ptsbma <= n*0.01 | n*0.99 <= ptsbma;
    %pp(:,7:9) = ptsbma(:,7:9) <= n*0.01; %for GRID    
    sum(pp(:,show))/length(pairs)  %COMPARE TO SHUFF
    pptsbma = pp;
    clear pp; clear n; clear show; 
%% 
ccofts = [ccof(ctsbma(:,1),ctsbma(:,2)), ccof(ctb, ctm),... dfdfsdf
ccof(ctsbma(:,4),ctsbma(:,5)), ccof(csb,csm)]
%atime= [mean(abs(ctsbma(:,1))), mean(abs(ctsbma(:,2))),mean(abs(ctb)),mean(abs(ctm)) ]
%aspace= [mean(abs(ctsbma(:,4))),mean(abs(ctsbma(:,5))), mean(abs(csb)), mean(abs(csm))]
    
%%    
%    ctsbma(43,:)=[]; pptsbma(43,:)=[]; pairs(43,:)=[]; ptsbma(43,:)=[];
    
%% COMPARE TO OUT GROUP
disp 'out group'
    tt = {ptcb,ptcm,pscb,pscm};
    ttt = {ctb, ctm, csb, csm};
    zptsbma = zeros(len(pairs),len(tt));
    for ii = 1: len(tt)
        %ctsbma(:,ii) = tt{ii}(:,1);
        %tt{ii}(:,2:end)=ttt{ii};
        for i = 1:len(pairs)
            t = [tt{ii}(i,1) ttt{ii}];
            %[~,I]=sort(abs(t),'descend');if ii>6;[~,I]=sort(t,'descend');end
            [~,I]=sort(t,'descend');
            pp = find(I==1); %index of non shuffled
            zptsbma(i,ii) = pp;
        end
    end
    n = len(ctb) + 1;
    pp = zptsbma<= n*0.01;  sum(pp(:,1:4))/length(zptsbma) %LARGEST VAL (FIRST DESC)
    pp = n*0.99 <= zptsbma; sum(pp(:,1:4))/length(zptsbma) %SMALLEST VAL (LAST DESC)
    pp = zptsbma <= n*0.01 | n*0.99 <= zptsbma; zpptsbma = pp; %zpptsbma(43,:)=[];
    sum(pp(:,:))/length(zptsbma) %COMPARE TO NON GROUP
    %sum(pptsbma(:,[1,2,4,5]))/len(pptsbma(:,1)) 


    
    
% %%
% ccofspace= [ccof(ctsbma(:,4),ctsbma(:,5)), ccof(zb,zm)]
% aspace= [mean(abs(ctsbma(:,4))), mean(abs(ctsbma(:,5))),mean(abs(zb)),mean(abs(zm)) ]

    
    
    
    %%
    par.n=1000;
    pp = ptsbma;
    pp(pp==1)=par.n*0.01;%for runs smaller than 100
    pp(pp>par.n*0.01)=0; %indices greater than 1
    pp(pp>0)=1;
    pptsbma = logical(pp);
    sum(ptsbma(:,:)<=par.n*0.01|ptsbma(:,:)==1)/length(pairs)*100    
    sum(pptsbma)/length(pairs)*100
    


    
    
    
    

a=0; aa=0;
for i = 1:len(cels)
    a(i)= cells{cels(i)}.before.gridscore;%/len(cels);
    aa(i)=cells{cels(i)}.midall.gridscore;%/len(cels);
    %if t==-2/len(cels); t=0; end;aa = aa +t;
end
stats(a) 
aa(aa==-2)=0;stats(aa)

a=0; aa=0;
for i = 1:len(cells)
    a(i)= cells{(i)}.before.gridscore;%/len(cels);
    aa(i)=cells{(i)}.midall.gridscore;%/len(cels);
    %if t==-2/len(cels); t=0; end;aa = aa +t;
end
a(a==-2)=0;stats(a) 
aa(aa==-2)=0;stats(aa)

%[133] bad cell, why? insignificant