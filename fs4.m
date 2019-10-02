%Sig non sig 
% f3s3
function s4(ctsbma,pptsbma)
cts=ctsbma;ppts=pptsbma;mfs=20;lfs=3; m='.';yls='space';
xl=mnmx(cts(:,[1 2])); yl=mnmx(cts(:,[4 5])); mss='markersize';
fig=figure(1002); clf; set(fig,'color','w', 'Position', [200 70 400 800]);
a=1; b=4; tt = 'pre'; subplot(211); f();
a=2; b=5; tt = 'dur'; subplot(212); f();

function f()
t= ppts(:,a)& ppts(:,b);plot(cts(t,a),cts(t,b),m,mss,mfs-3.2,'linewidth',lfs-1.5);tss=sum(t); 
hold on;
t= ppts(:,a)&~ppts(:,b);plot(cts(t,a),cts(t,b),m,mss,mfs-3.4,'linewidth',lfs-1.5);tsn=sum(t);
t=~ppts(:,a)& ppts(:,b);plot(cts(t,a),cts(t,b),m,mss,mfs-3.6,'linewidth',lfs-1.5);tns=sum(t);
t=~ppts(:,a)&~ppts(:,b);plot(cts(t,a),cts(t,b),m,mss,mfs-3.8,'linewidth',lfs-1.5);tnn=sum(t);
title(sprintf('%s: corrs time vs %s shuff',tt,yls)); [~,~,s]=ccof(cts(:,a),cts(:,b),2);
text(0.05,0.9,s,'Units','normalized');
axis square; hold off; xlabel('time'); ylabel(yls); axis 'tight';xlim(xl);ylim(yl);
legend({['sig sig n=' n2(tss)]; ['sig non n=' n2(tsn)]; ...
        ['non sig n=' n2(tns)]; ['non non n=' n2(tnn)]},'location','southeast')
end
end

% plot(cts( ppts(:,a)& ppts(:,b),a), cts( ppts(:,a)& ppts(:,b),b),m,mss,mfs-3.2,'linewidth',lfs-1.5); 
% hold on;
% plot(cts( ppts(:,a)&~ppts(:,b),a), cts( ppts(:,a)&~ppts(:,b),b),m,mss,mfs-3.4,'linewidth',lfs-1.5);
% plot(cts(~ppts(:,a)& ppts(:,b),a), cts(~ppts(:,a)& ppts(:,b),b),m,mss,mfs-3.6,'linewidth',lfs-1.5);
% plot(cts(~ppts(:,a)&~ppts(:,b),a), cts(~ppts(:,a)&~ppts(:,b),b),m,mss,mfs-3.8,'linewidth',lfs-1.5);
% a=1; b=4; %a=2; b=5;
% 'sig:sig   sig:non   non:sig   non:non'
% [sum(ppts(:,a)&ppts(:,b)),...
% sum(ppts(:,a)&~ppts(:,b)),...
% sum(~ppts(:,a)&ppts(:,b)),...
% sum(~ppts(:,a)&~ppts(:,b)) ]
% if(a==1) disp('pre'); else disp('dur'); end

%Pre=  bar([40 27 17 15; 1 37 0 61]');
%Dur=  bar([1, 37, 0, 61])
% title('significance [temp:xgs]')
%xlabel(['sig:sig   ', 'sig:non   ', 'non:sig   ', 'non:non']);
%legend({'PRE','DUR'})