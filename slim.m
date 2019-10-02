%set x an y limits equal based on respective limits to encompass all data
function [ lim ] = slim( axa )
mnx = [];
for i = 1:len(axa)
    mnx(end+1,:) = axa(i).XLim;
    mnx(end+1,:) = axa(i).YLim;
end
lim = minmax(mnx(:)');
end


