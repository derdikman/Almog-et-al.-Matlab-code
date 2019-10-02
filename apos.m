function p = apos(ncol, naxes,ai)
    nRow = ceil( naxes / ncol ) ;
    rowH = 0.70 / nRow ;  colW = 0.70 / ncol ; %make smaller for more label room
           %offset
    colX = 0.05 + linspace( 0, 0.90, ncol+1 ) ; colX = colX(1:end-1) ; %.05 .90
    rowY = 0.05 + linspace( 0.90, 0, nRow+1 ) ; rowY = rowY(2:end);
    rowId = ceil( ai / ncol );
    colId = ai - (rowId - 1) * ncol;
    x = colX(colId); y = rowY(rowId);
    p = [x, y, colW, rowH];
end

