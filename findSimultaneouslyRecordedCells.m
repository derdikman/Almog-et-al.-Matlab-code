function [groups, cells] = find_simultaneously_recorded_cells(cells)
    groups = [];
    for i=1:length(cells)
        if ~isfield(cells{i}, 'midall')
            cells{i}.midall = midAll(cells{i});
        end
        a = cells{i};
        if ~isempty(a.midall) %skip cells with non middle %middle
            key = sprintf('g_%s_%s', a.id, a.date);
            if isfield(groups, key)
                t = groups.(key);
                t = [t a];
                groups.(key) = t;
            else
                groups = setfield(groups, key, a);
            end
        end
        %    fprintf('%d: bad cell %f rmthr %f', i, a.max_r, rate_map_thresh);
    end
    
    %validate
    keys = fieldnames(groups);
    easy = {};
    for i = 1:length(keys)
        key = keys{i};
        group = groups.(key);
        easy{i} = group;
        t = group(1);
        t = length(t.before.px);
        for j = 1:length(group) %CHECK TIMES !!
            tt = group(j); %cells{group(j)};
            tt = length(tt.before.px);
            if t ~= tt
                fprintf('%f %f whats up with this group b.length %s?\n',t, tt, key);
            end
        end
    end
    groups = easy;
    
    for j = 1:length(groups)
        g = groups{j}; 
        t = [];
        for i = 1:length(g)
            
        end
    end
    
    
    %     for j = 1:length(groups)
    %         g = groups{j};
    %         t = [];
    %         for i = 1:length(g)
    %             c = g(i);
    %             for mid = 1:length(c.middle);
    %                 t = [t c.middle{mid}.pt(1)];
    %             end
    %         end
    %         t = unique(t);
    %         t = sort(t)';
    %         assert(min(diff(t))/60>=1,'bins not even 1min apart');
    %         for i = 1:length(g)
    %             c = g(i);
    %             for mid = 1:length(c.middle);
    %                 %t = [t c.middle{mid}.pt(1)];
    %             end
    %         end
    %        % groups{j} = g;
    %     end
    
    %ADD CODE TO ALIGN BINS
    
    %fprintf('grouped ADD CODE TO ALIGN BINS\n');
end

function midall = midAll(c)
    
    
    midall.px = []; midall.py = []; midall.pt = [];
    midall.sx = []; midall.sy = []; midall.st = [];
    for mid = 1:length(c.middle)
        midall.px = [midall.px c.middle{mid}.px];
        midall.py = [midall.py c.middle{mid}.py];
        midall.pt = [midall.pt c.middle{mid}.pt];
        midall.sx = [midall.sx c.middle{mid}.sx];
        midall.sy = [midall.sy c.middle{mid}.sy];
        midall.st = [midall.st c.middle{mid}.st];
    end
    midall.px = double(midall.px);midall.py = double(midall.py);midall.pt = double(midall.pt);
    midall.sx = double(midall.sx);midall.sy = double(midall.sy);midall.st = double(midall.st);
    midall.rm =  Create_Rate_Map(midall.px, midall.py, midall.pt,...
        midall.sx, midall.sy, midall.st);
    midall.ac = Cross_Correlation(midall.rm, midall.rm);
    midall.max_r = max(midall.rm(:));
    midall.gridscore = gridscore(midall.ac, c.ind);
    midall.module = Find_Module(midall.ac);
    midall.exists = false;
    if midall.max_r ~= 0 && midall.max_r ~= 50 && midall.gridscore ~= -2
        midall.exists = true;
    end
end
