function o = plotARP(x,y,arg)

rd=2;ax=gca;show=[1 1 1];tx=0.1; ty=0.9;lw=2; ms=10; 
slmd=0; setlim=true;


if isfield(arg,'show');  show=arg.show; end %[s s s] which values to show
if isfield(arg,'setlim');  setlim=arg.setlim; end
if isfield(arg,'slmd');  slmd=arg.slimd; end %sets limits by data
if isfield(arg,'rd');    rd=arg.rd; end %round by
if isfield(arg,'ax');    ax=arg.ax; end %axes
if isfield(arg,'tx');    tx=arg.tx; end %location of text
if isfield(arg,'ty');    ty=arg.ty; end
if isfield(arg,'lw');    lw=arg.lw; end %line width
if isfield(arg,'ms');    ms=arg.ms; end %markser size
%if isfield(arg,'');   =arg.; end

o=plot(ax,x,y,'.','markersize',ms); hold on;

l=slim(ax); if slmd; l=slimd(ax); end 
if setlim
ax.XLim = l; ax.YLim = l;
end
axis(ax,'square');

f = fit(x,y,'poly1');
plot(ax.XLim,f(l),'--','linewidth',lw); 

a=f.p1; [r, p]=ccof(x,y); as=''; rs=''; ps=''; 
ff=['%0.' n2(rd) 'f'];
if show(1); as= sprintf(['a=' ff], round(a,rd)); end
if show(2); rs= sprintf(['r=' ff], round(r,rd)); end
%if show(3); ps= sprintf(['p=' ff], round(p,rd)); end
if show(3); ps= pstr(p,rd); end

text(tx,ty,[as ' ' rs ' ' ps],'Units','normalized');

end