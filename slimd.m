%set x an y limits equal based on respective data to encompass all data
function [ lim ] = slimd( axa )
mnx = [];
for i = 1:len(axa)
%     ax = axa(i);
    for ii=1:len(axa(i).Children)
        mnx(end+1,:) = minmax(axa(i).Children(ii).XData);
        mnx(end+1,:) = minmax(axa(i).Children(ii).YData);
    end
    %d = get(ax.Children');
%     mnx(end+1) = min(ax.Children(2).XData);
%     mnx(end+1) = max(ax.Children(2).XData);
%     mnx(end+1) = min(ax.Children(2).YData);
%     mnx(end+1) = max(ax.Children(2).YData);
end
% lim(1) = min(mnx);
% lim(2) = max(mnx);
lim=minmax(mnx(:)');
end


