function a = patchTrajectoryLinear(pt,px,py,dt,dtthresh)
    pt = round(pt,3); dt = round(dt,3); dtthresh = round(dtthresh,3);
    if isnan(dt)
            a.x = 0; 
            a.y = 0; 
            a.t = 1;
        return
    end
    t = zeros(ceil(pt(end)/dt),1); x = t; y = t; ttt=[];
    t(1) = pt(1); x(1) = px(1); y(1) = py(1);
    ii = 1;
    for i = 2:length(pt)
       %if time from previous step is larger than thresh
       if pt(i) - pt(i-1) > dtthresh %%pt(i) - pt(i-1)
            it = i;
            tt = pt(i-1):dt:pt(i); %will always include first number, then all numbers counting 
                                   %forward that are less than or equal to limit number             
            tt = tt(2:end-1);%use only times in between limits %WAS end-1)
            step = dt/(pt(i)- pt(i-1)); %number < 1, shouldn't be based on thresh?
            assert(abs(step)<1);
            r = (1:length(tt))*step;
            tx = px(i-1) + r*(px(i)- px(i-1));
            ty = py(i-1) + r*(py(i)- py(i-1));
            
            %gap not large enough to add timestep
            if isempty(tt)
                tt = pt(i); 
                tx = px(i);
                ty = py(i);
            %make sure always use actual data
            elseif tt(end)==pt(i)
                tx(end) = px(i);
                ty(end) = py(i);
            else
                tt(end+1) = pt(i);
                tx(end+1) = px(i);
                ty(end+1) = py(i);
            end
            %remove mini timesteps
            if len(tt)>1 && tt(end)-tt(end-1)<dt
                tt(end-1)=[];tx(end-1)=[];ty(end-1)=[];
                
            end
                
            
       else
            tt = pt(i);        
            tx = px(i);
            ty = py(i);
       end
       ii = ii+length(tt);
       t(ii-length(tt)+1:ii) = tt;%t = [t tt];
       %ttt = [ttt tt];
       x(ii-length(tt)+1:ii) = tx;%x = [x tx];
       y(ii-length(tt)+1:ii) = ty;%y = [y ty];
       assert(t(ii)-t(ii-1)<2*dtthresh)
    end
    a.x = x(1:ii); 
    a.y = y(1:ii); 
    a.t = t(1:ii);
    assert(min(diff(a.t))>dt/2);
    assert(max(diff(a.t))<2*dtthresh,'max(diff(t))<dtthresh');
    assert(sum(a.t==0)<=1,'sum(isnan(a.x))==0 trajectory');
    assert(len(a.t)==len(a.x)&&len(a.x)==len(a.y));
    
%     is = 1:length(pt);
%     t = diff(pt);
%     is = is(t>dtthresh);
%     patches = cell(1,length(is));
%     for i = 1:length(is)
%         %i = 441
%         t = pt(is(i)):dt:pt(is(i)+1);
%         patches(i).t = t(2:end-1);
%         step = dt/(pt(is(i)+1)- pt(is(i))) %number < 1
%         assert(abs(step)<1);
%         r = (1:length(t))*step;
%         patches(i).x = px(is(i)) + r*(px(is(i)+1)- px(is(i)));
%         patches(i).y = py(is(i)) + r*(py(is(i)+1)- py(is(i)));
%         [ c.px(is(i)) x c.px(is(i)+1)]
%         [ c.py(is(i)) y c.py(is(i)+1)]
%         [ c.pt(is(i)) t c.pt(is(i)+1)]
%         diff([ c.px(is(i)) x c.px(is(i)+1)])
%         diff([ c.py(is(i)) y c.py(is(i)+1)])
%         diff([ c.pt(is(i)) t c.pt(is(i)+1)])
    
end

