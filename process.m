function r = process(data, binsizeminutes)
nbins = 50;
r.ind = data.ind;
ep.x1 = 0; ep.x2 = 0; ep.y1 = 0; ep.y2 = 0; ep.t = 1;
es.x = 0; es.x2 = 0;  es.y = 0;  es.y2 = 0; es.ts =1;
%before makeSession=createSession?
r.before = createSession(data.B(1).pos_data, data.B(1).spike_data,nbins,[-inf inf]);
%middle
r.midall = createSession(data.B(2).pos_data, data.B(2).spike_data,nbins,[0,binsizeminutes*60]);
%after
if length(data.B) == 3 
    r.after = createSession(data.B(3).pos_data, data.B(3).spike_data,nbins,[-inf inf]);
else
    r.after = createSession(ep, es, nbins,[-inf inf]);
end
r.tet  = data.tetrode;
r.cell  = data.cell;
r.type = data.cell_type_obj_new;
r.area = data.area;
r.id   = data.rat;
r.date = data.date;
end

% bins = bin_trajactory(tr, binLength*60, maxTime*60); %BIN LENGTH
% for i = 1:length(bins)
%     tr = bins{i};
%     [bins{i}.rm, bins{i}.max_r] = Create_Rate_Map(tr.px, tr.py, tr.pt, tr.sx, tr.sy, tr.st,true,nbins);
%     bins{i}.ac = Cross_Correlation(bins{i}.rm, bins{i}.rm);
%     bins{i}.ac2 = xcorr2(bins{i}.rm); 
%     %fprintf('%d %d\n', r.ind, i);
%     %bins{i}.gridscore = gridscore(bins{i}.ac, r.ind);
%     bins{i}.gridscore = gridscore2(bins{i}.ac2, r.ind);
%     bins{i}.module = Find_Module(imgaussfilt(bins{i}.ac2, 2,'FilterDomain','spatial'));
%     bins{i}.exists = false;
%     if bins{i}.max_r ~= 0 %&&  bins{i}.gridscore ~= -2 %bins{i}.max_r ~= 50 &&
%         bins{i}.exists = true;
%     end
% end
% %r.middle = bins; %CHANGE BACK
% r.midall = bins{1};