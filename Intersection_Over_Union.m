function iu=Intersection_Over_Union(x1,y1,x2,y2)   
    if isreal(x1) && isreal(x2)&& isreal(y1) && isreal(y2)
        [xu, yu] = polybool('union', x1, y1, x2, y2);
        [xi, yi] = polybool('intersection', x1, y1, x2, y2);
        intersection = polyarea(xi,yi);
        union = polyarea(xu,yu);
    else
        intersection = 1;
        union=-1;
    end
    iu=intersection/union;
    %figure();   clf;title(iu); hold on; plot(xa,ya);plot(xb,yb); %plot(x1,y1);plot(x2,y2);
end