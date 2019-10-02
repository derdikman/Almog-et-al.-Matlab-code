function o = imgsc(a,sig,ydirfalse)

    if nargin == 1
        o=imagesc(a); set(gca,'ydir','norm');
    elseif nargin ==2
        o=imagesc(imgaussfilt(a,sig,'FilterDomain','spatial')); set(gca,'ydir','norm');

    else 
        if sig==0
            o=imagesc(a);
        else
            o=imagesc(imgaussfilt(a,sig,'FilterDomain','spatial'));
        end
    end
    
    axis(gca,'square'); colormap jet;


end