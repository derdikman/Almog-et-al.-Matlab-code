%c session
%str strings for labels
%ax axis
function ax = plotRM(ax,c,str)%m,n,l    
    title(ax, str.t);xlabel(ax, str.x);ylabel(ax, str.y);hold(ax, 'on');
    
    imagesc(ax, imgaussfilt(c.rm,2,'FilterDomain','spatial')); 
    
    axis(ax,'tight'); set(ax,'ydir','normal');colormap(ax,'jet'); axis(ax,'square');
end

