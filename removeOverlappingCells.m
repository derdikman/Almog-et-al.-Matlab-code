function [gids, cids] = removeOverlappingCells(groups, gids, cids, threshRemoveSpikes, threshRemoveCell, blacklist)
    rgs = [];
    percentremoved=[];
    for ri = 1:length(gids)
        gis = gids(ri); cis = cids{gids(ri)};%gis}
        g = groups{gis};%gis};
        g = g(cis);%cis);
        rcs = [];
        for j = 1:length(g)-1
            for k = j+1:length(g)
                c1 = g(j).before; c2 = g(k).before;
                s1=floor(c1.st*1000)+1; s2=floor(c2.st*1000)+1; min_time=min(min(s1),min(s2)); %TIME IN MS
                offset = 0; s1 = s1-min_time+offset; s2 = s2-min_time+offset; max_time=max(max(s1),max(s2))+offset;
                [s1i, s2i] = removeOverlappingSpikes(s1,s2, 1); rs1=s1(s1i); rs2=s2(s2i);
                %fprintf('####b.g%2d i%3dxi%3d: %3d :: %3d\n',ri, g(j).ind ,g(k).ind,length(s1)-length(rs1),length(s2)-length(rs2));
                rem = length(s1)-length(rs1); %assert(abs( (length(s1)-length(rs1)) - (length(s2)-length(rs2)) )<3)
                brem = max(rem/length(s1),rem/length(s2));
                if brem>threshRemoveCell || false
                     fprintf('b.g%2d: spikes removed %3dx%3d = %3d, --- c1 %2d%% c2 %2d%% \n',gids(ri), g(j).ind,g(k).ind,rem,...
                     100*round(1-length(rs1)/length(s1),2),100*round(1-length(rs2)/length(s2),2));
                else
                    percentremoved(end+1)=brem;
                end
                if brem >  threshRemoveCell;
                    %remove cell with higher percentage spikes removed
                    rid = cis(j);
                    %s2 got reduced more
                    if length(rs2)/length(s2) < length(rs1)/length(s1)  
                        rid = cis(k);
                    end
%                     Sfprintf('g.%2d: b c.%3d: %3f%%\n',gids(ri), rid, 100*max(brem,mrem));
                    rcs = [rcs rid];
                end
                c1 = g(j).midall; c2 = g(k).midall;
                s1=floor(c1.st*1000)+1; s2=floor(c2.st*1000)+1; min_time=min(min(s1),min(s2)); %TIME IN MS
                offset = 0; s1 = s1-min_time+offset; s2 = s2-min_time+offset; max_time=max(max(s1),max(s2))+offset;
                [s1i, s2i] = removeOverlappingSpikes(s1,s2, 1); rs1=s1(s1i); rs2=s2(s2i);
                rem = length(s1)-length(rs1); %assert(abs( (length(s1)-length(rs1)) - (length(s2)-length(rs2)) )<3)
                mrem = max(rem/length(s1),rem/length(s2));
                if mrem>threshRemoveCell || false
                     fprintf('m.g%2d: spikes removed %3dx%3d = %3d, --- c1 %2d%% c2 %2d%% \n',gids(ri), g(j).ind,g(k).ind,rem,...
                     100*round(1-length(rs1)/length(s1),2),100*round(1-length(rs2)/length(s2),2));
                else
                    percentremoved(end+1)=mrem;
                end
                if mrem > threshRemoveCell
                    %remove cell with higher percentage spikes removed
                    rid = cis(j);
                    %s2 got reduced more
                    if length(rs2)/length(s2) < length(rs1)/length(s1)  
                        rid = cis(k);
                    end
                    %fprintf('g.%2d: m c.%3d: %3f%%\n',gids(ri), rid, 100*max(brem,mrem));
                    rcs = [rcs rid];
                end
                %hacky...
                if ismember(g(j).ind,blacklist);rcs = [rcs cis(j)];end
                if ismember(g(k).ind,blacklist);rcs = [rcs cis(k)];end
                %disp('__');
            end
        end
        %remove cells
        if ~isempty(rcs)
            rcs = unique(rcs);
            fprintf('$ removing %d cells from g.%d\n',length(rcs),gids(ri));
            %rcs
            cids{gids(ri)} = setdiff(cids{gids(ri)},rcs);
            
            if len(cids{gids(ri)}) < 2 %only groups with at least 2
                'here'
                rgs = [rgs ri];
            end
        end
    end
    %%
    fprintf('removing %d groups\n',length(rgs));
    gids(rgs)=[];
    
end