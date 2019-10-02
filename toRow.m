function x = toRow(x)
    if iscolumn(x)
        x = x';
    end
end