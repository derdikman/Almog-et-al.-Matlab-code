%set x an y limits equal based on respective data to encompass all data
function [ lim ] = slimdy( axa )
mnx = [];
for i = 1:len(axa)
    for ii=1:len(axa(i).Children)
        mnx(end+1,:) = minmax(axa(i).Children(ii).YData);
    end
end
lim=minmax(mnx(:)');
end
