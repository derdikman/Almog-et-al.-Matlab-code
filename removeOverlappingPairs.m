 function pairsr = removeOverlappingPairs(cells,pairst,threshRemoveSpikes, threshRemoveCell)

percentremoved=[];
rcs = [];
for i = 1:len(pairst)
    cc1 = cells{pairst(i,1)}; cc2=cells{pairst(i,2)};
    c1 = cc1.before; c2 = cc2.before;
    s1=floor(c1.st*1000)+1; s2=floor(c2.st*1000)+1; 
    [s1i, s2i] = removeOverlappingSpikes(s1,s2, threshRemoveSpikes); rs1=s1(s1i); rs2=s2(s2i);
    fprintf('i%3dxi%3d: %3d :: %3d\n', cc1.ind ,cc2.ind,len(s1)-len(rs1),len(s2)-len(rs2));
    rem = len(s1)-len(rs1); %assert(abs( (len(s1)-len(rs1)) - (len(s2)-len(rs2)) )<3)
    brem = max(rem/len(s1),rem/len(s2));
    if brem > threshRemoveCell
        fprintf('bspikes removed %3dx%3d = %3d, --- c1 %2d%% c2 %2d%% \n', cc1.ind,cc2.ind,rem,...
            100*round(1-len(rs1)/len(s1),2),100*round(1-len(rs2)/len(s2),2));
        rcs = [rcs i];
        
    else %NOT REMOVED BY BEFORE
        c1 = cc1.midall; c2 = cc2.midall;
        s1=floor(c1.st*1000)+1; s2=floor(c2.st*1000)+1;
        [s1i, s2i] = removeOverlappingSpikes(s1,s2, threshRemoveSpikes); rs1=s1(s1i); rs2=s2(s2i);
        rem = len(s1)-len(rs1); %assert(abs( (len(s1)-len(rs1)) - (len(s2)-len(rs2)) )<3)
        mrem = max(rem/len(s1),rem/len(s2));
        if mrem>threshRemoveCell || false
            fprintf('mspikes removed %3dx%3d = %3d, --- c1 %2d%% c2 %2d%% \n', cc1.ind,cc2.ind,rem,...
                100*round(1-len(rs1)/len(s1),2),100*round(1-len(rs2)/len(s2),2));
             rcs = [rcs i];
        end
        
    end
end

pairst(rcs,:) = [];
pairsr = pairst;

end
