warning off;
%%%
binsize = 45;
%fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_c_midall_gridscore.mat',binsize); first set
fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_d_patchtraj_rayleigh',45);
load(fn);
[groups ~] = findSimultaneouslyRecordedCells(cells);
gridThreshBef = 0.3; gridThreshMid = 0.25;
[bgids, bcids] = groupsByGridThresh(groups, gridThreshBef, gridThreshMid);
gids = []; cids = [];
[gids, cids] = removeOverlappingCells(groups, bgids, bcids, -1, 0.07,[]); %why 133 bad?
n = 250;
%grid 
gsshufb = [];
gsshufm = [];
gsshufa = [];
ii = 0;
tic
%shuffles gridscore
for ri = 1:length(gids)
    gis = gids(ri); cis = cids{gids(ri)};%gis}
    g = groups{gis};%gis};
    g = g(cis);%cis);
    for ci = 1:length(g)
        ii = ii+1; 'ind'; g(ci).ind
        %%%%%% midall %%%%%%
        c = g(ci).before;
        maxt = max(c.pt); disc = mean([c.pt(1:end-1) c.pt(2:end)],2);
        d = maxt/n; t = zeros(1,n);
        for hi = 1:n %n+1 because wi = (i-1)*d; to include 0 shift
            %fprintf('done %.2f%%\n',round(hi/(n+1),2));
            toff = (hi-1)*d; st = c.st+toff; st(st>maxt) = st(st>maxt)-maxt;
            assert(max(st)<maxt);
            si = discretize(st,[-Inf;disc; +Inf]);sx = c.px(si); sy = c.py(si);
            [rateMap, ~] = createRateMap(c.px, c.py, c.pt, sx, sy, st, true, 50);
            ac = xcorr2(rateMap); t(hi) = gridscore2(ac, 2);
        end% [m,mi] = max(t)
        gsshufb(ii,:) = [g(ci).ind,t];
        %%%%%% midall %%%%%%
        c = g(ci).midall;
        maxt = max(c.pt); disc = mean([c.pt(1:end-1) c.pt(2:end)],2);
        d = maxt/n; t = zeros(1,n);
        for hi = 1:n %n+1 because wi = (i-1)*d; to include 0 shift
            %fprintf('done %.2f%%\n',round(hi/(n+1),2));
            toff = (hi-1)*d; st = c.st+toff; st(st>maxt) = st(st>maxt)-maxt;
            assert(max(st)<maxt);
            si = discretize(st,[-Inf;disc; +Inf]);sx = c.px(si); sy = c.py(si);
            [rateMap, ~] = createRateMap(c.px, c.py, c.pt, sx, sy, st, true, 50);
            ac = xcorr2(rateMap); t(hi) = gridscore2(ac, 2);
        end% [m,mi] = max(t)
        gsshufm(ii,:) = [g(ci).ind,t];
        %%%%%% after %%%%%%
        c = g(ci).after; t = zeros(1,n);
        if len(c.st) > 1
            maxt = max(c.pt); disc = mean([c.pt(1:end-1) c.pt(2:end)],2);
            d = maxt/n; 
            for hi = 1:n %n+1 because wi = (i-1)*d; to include 0 shift
                %fprintf('done %.2f%%\n',round(hi/(n+1),2));
                toff = (hi-1)*d; st = c.st+toff; st(st>maxt) = st(st>maxt)-maxt;
                assert(max(st)<maxt);
                si = discretize(st,[-Inf;disc; +Inf]);sx = c.px(si); sy = c.py(si);
                [rateMap, ~] = createRateMap(c.px, c.py, c.pt, sx, sy, st, true, 50);
                ac = xcorr2(rateMap); t(hi) = gridscore2(ac, 2);
            end% [m,mi] = max(t)
        end
        gsshufa(ii,:) = [g(ci).ind,t];
        %%
    end
end
toc

gvbp = [];
for i = 1: size(gsshufb,1);
t = gsshufb(i,:);
[~, mi] = max(t(2:end));
gvbp(i,:) = [t(1), mi];
%[s,I] = sort(abs(t),'descend');
%pb = round(find(I==1)/length(I),2); %index of non shuffled
end
unique(gvbp(:,2))

gvbm = [];
for i = 1: size(gsshufm,1);
    t = gsshufm(i,:);
    if t(2)-max(t(3:end)) > 0 %how much is gridscore smaller than shuffled
        t(1)
        disp(round(t(2)-max(t(3:end)),2))
    end
    [m, mi] = max(t(2:end));
    if m == -2
        mi = 100;
    end
    %[s,I] = sort(abs(t),'descend');
    %pb = round(find(I==1)/length(I),2); %index of non shuffled
    gvbm(i,:) = [t(1), mi];
end
unique(gvbm(:,2))


gvba = [];
for i = 1: size(gsshufb,1);
t = gsshufa(i,:);
[~, mi] = max(t(2:end));
gvba(i,:) = [t(1), mi];
%[s,I] = sort(abs(t),'descend');
%pb = round(find(I==1)/length(I),2); %index of non shuffled
end
unique(gvba(:,2))

what are shifting parameters max st or pt
