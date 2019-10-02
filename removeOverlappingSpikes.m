%takes as input two arrays of spike times in chronological order, overlap in same time units
%deletes BOTH spikes in case of overlap?
%returns INDICES of spikes to REMOVE
function [s1keep, s2keep] = removeOverlappingSpikes(s1,s2, overlap)

% %     c2 = cells{71}.before; c1 = cells{75}.before;
% %     s1=floor(c1.st*1000)+1; s2=floor(c2.st*1000)+1; mint=min(min(s1),min(s2)); %TIME IN MS
% %     s1 = s1-mint+1; s2 = s2-mint+1; maxt=max(max(s1),max(s2));
%     %for debugging
%     train1=zeros(1,maxt);train2 =zeros(1,maxt); train1(s1)=1; train2(s2)=1;%to compare after;
%     figure;plot(1:length(train1),[train1+1; train2]); %to compare after;
%    t1 = length(s1); t2 = length(s2); %to compare after; %debug
    s1orig = s1; s2orig = s2;
    s2i = 1; %index of current spike in s2
    s2Time = s2(1); %time of current spike in s2
    s2deletedTime = -inf;  %spike time of spike in s2 removed
    rem = 1; % count of spikes deleted
    for s1i = 1:length(s1) % index of current spike in s1
        s1Time = s1(s1i); %current spk1 time
        %this spike s2td in s2 already got deleted and was in overlap of previous s1 spike
        if s1Time <= s2deletedTime + overlap  %CHANGED FROM <
            %delete spikes overlapping with those already deleted in s2
            s1(s1i) = 0;
            %disp('here');
            rem = rem+1;%count spks deleted
            %this section is here because we are iterating by spk1 but the
            %last spike deleted in s2 could overlap with current spike and
            %so this spike should be deleted as well
        end
        %spk2 is within +- overlap
        %spk2time is smaller than overlap end
        while s2Time <= s1Time + overlap %CHANGED FROM <
            %spk2 greater than overlap beginning
            if s2Time >= s1Time - overlap %CHANGED FROM >
                %fprintf('%d %d\n',i,j);%spkt1(i),spkt2(j));
                s2deletedTime = s2Time; %last spike in s2 deleted;
                s1(s1i) = 0; s2(s2i) = 0; %deletes spike in spkt1 and 2 in overlap
                rem = rem+1;
            end
            s2i = s2i+1; %next spk2
            if s2i <= length(s2)
                s2Time = s2(s2i);
            else %only when finished with spk2 array
                s2Time = inf;
            end
        end
    end
   
    s1keep = s1~=0; s2keep = s2~=0;
    s1 = s1(s1~=0); s2 = s2(s2~=0);
% %     %     %for debugging
%      fprintf('spikes removed %dx%d = %3d, c1 %2d%% c2 %2d%% \n',1,2,rem,...
%          100*round(1-length(s1)/length(s1orig),2),100*round(1-length(s2)/length(s2orig),2));

%     train1=zeros(1,maxt);train2 =zeros(1,maxt); train1(s1)=1; train2(s2)=1;
%     figure;plot(1:length(train1),[train1+1; train2]);
end

