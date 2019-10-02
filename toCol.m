function x = toCol(x)
    if ~iscolumn(x)
        x = x';
    end
end

