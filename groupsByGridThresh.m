function [gids, cids] = groupsByGridThresh(groups, gridThreshBef, gridThreshMid)
    gids = []; cids = {};
    for i = 1:length(groups)
        g = groups{i};
        t = [];
        for j = 1:length(g)
            if g(j).before.gridscore > gridThreshBef...
            && g(j).midall.gridscore < gridThreshMid
            %if g(j).before.gridscore < g(j).midall.gridscore 
                t(end+1) = j;
            end
        end
        if length(t) >= 2 %only groups with at least 2
            cids{i} = t; 
            gids(end+1) = i;
        end
    end
end