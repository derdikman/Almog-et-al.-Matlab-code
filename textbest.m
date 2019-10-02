function hOut = textbest(ax,str,lvar,tvar)

l = legend(ax,str,lvar{:});
x = l.Position(1);y = l.Position(2);
delete(l);
text(ax,x,y,str,'Units','normalized',tvar{:});
if nargout
    hOut = t;
end
end