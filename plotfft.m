function f = plotfft(ax,x,fs,fmax,fdel,win) %fdel was 5 r 1000
    y=abs(fft(x-mean(x))); y = y(1:floor(length(x/2)));
    if exist('r','var')
        y = decimate(y,r); %y=abs(y); r = 1000;
    else r = 1; end;
    l = length(y);
    f = fs*(1:len(y))/l; %/2??
    if exist('fmax','var')
        f = f(f<=fmax+1);
    end
    %deletes first few; median 1000
    y = y(1:len(f));
    y(f<fdel)=mean(y(f>=fdel));
    %y = movmean(y,10);
    if exist('win','var') && round(win*fs/1000) > 0
        y(f>=fdel) = movmean(y(f>=fdel),round(len(y(f>=fdel))*win/(r*1)));
    end
    y = y-min(y); y = y/max(y);
    %p = plot(ax,f,abs(y(1:floor(l/2))));
    p = plot(ax,f,y,'linewidth',2); ylabel(ax,len(y));
    %plot(ax,f,[toCol(real(y(1:floor(l/2)))) , toCol(complex( y(1:floor(l/2) ) ) ) ] );
    %hold(ax,'on'); %plot(ax,[7.9999 8],[ylim(ax)],'r--','linewidth',1.2); hold(ax,'off'); %max(y(1:round(l/2)) )
    ylim(ax,[0 1]);
    xlim(ax,[0 max(f)]);
    xlabel(ax,'Hz');
end

%{

  y=fft(train-mean(ftrain));y(1:500)=0; %deletes first few;
         l = length(y); l2 = floor(l/2);
         y = abs(y(1:l2));
         f = 1000*(1:length(y))/(length(y)*2);
        % y = smooth(y,1e2);%'sgolay',
         fi = (length(y)*    20    *2/1000);
         plot(ax,f(1:fi),smooth(y(1:fi),500),'linewidth',3); ax.YLim = [0 100]; hold(ax,'on'); 
         plot(ax,[7.9999 8],[ylim(ax)],'r.','linewidth',1.2,'markersize',15);  %max(y(1:round(l/2)) )
         xlabel(ax,'Hz');
         title(ax,sprintf('gs %.1f  len %d ',round(c.gridscore,1),length(y)),'color','m','fontsize',tfs); ax.YLim = [0 100]; %hold(ax,'off'); 

%}