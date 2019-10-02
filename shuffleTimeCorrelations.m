function pcor0 = shuffleTimeCorrelations (c1, c2, p)
    [t1, t2] = createMsSpikeTrain(c1.st, 100, c2.st); %good, leading 0's lower numbers
    %t1=createMsSpikeTrain(c1.st);
    %t2=createMsSpikeTrain(c2.st);
    %delta to shift train by
    d = floor(length(t1)/p.n) ;%right?  %amount to shift by each step
    pcor0 = zeros(1,p.n);
    i=0;
    while i < p.n %n+1 because wi = (i-1)*d; to include 0 shift
        i=i+1;
        %sprintf('done %.2f\%',round(i/(p.n+1),2));
        wi = (i-1)*d; %first one calculated is 0
        shiftedi = [wi+1:length(t1) 1:wi]; %cyclic
        t2shifted = t2(shiftedi);
        pcor = timeCorrelationSmoothed(t1,t2shifted,p);
        pcor0(i) = pcor;
        if sum( abs(pcor0(2:end)) >= abs(pcor0(1)) ) >= ceil(p.n*p.pval)
                %['time abs(corr(1)) not significant shuffling i=' n2(i) ' pval=' n2(p.pval) ' more than ' n2(ceil(p.n*p.pval)) ' or gs(1)= '  n2(pcor0(1))]
           i=inf;
        elseif i==p.n
            ['TIME SIGNIFICANT' n2(pcor0(1))]
        end
    end
end

    %disp(pcor0);
    %[mx, imx] = max(abs(pcor0))
    %time_scale=( (1:length(p{i,1})) - ((length(p{i,1})-1)/2) - 1)/1000; %goes from -lag 0 +lag in ms
    %figure;
     %plotmatrix(time_scale',cell2mat(p(:,1))','-');
     %figure
     %plotmatrix((1:maxt)',trains','-');
    %plotmatrix(cell2mat(p(1,2))',cell2mat(p(:,1))')