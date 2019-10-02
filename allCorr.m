function A = allCorr(a, b) %
    A = zeros(length(a));
    for i = 1:length(a)
        for j = 1:length(b)
            t1 = a{i}(:); t2 = b{j}(:);
            %all this work to make two non-nan same length vectors!
            t1(isnan(t1))=0; t2(isnan(t2))=0;
            if length(t1) ~= length(t2)
                if mod(length(t1) + length(t2), 2) ~= 0
                    t2 = [t2 ; 0];
                end
                t1 = padarray(t1, [ceil((max([length(t1) length(t2)]) - (length(t1)))/2) 0]);
                t2 = padarray(t2, [ceil((max([length(t1) length(t2)]) - (length(t2)))/2) 0]);
            end
            A(i,j) = corr(t1, t2); 
        end
    end
end