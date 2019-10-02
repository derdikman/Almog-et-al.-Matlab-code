a = zeros(1000,1); l = length(a);b = a; clf;
%d = 3;
%b(1:d:length(a)) = 1;
e= 100.5; 
b(e:e:length(a)) = 1;

f= fft(b);%fft(b); %/length(a);   fwht
figure(4);
subplot(311);
%plot(f);
plot(abs(f(1:round(length(a)/2)))); title(num2str(round(length(a)/e))); hold on; plot(round(length(a)/e),0,'r*'); hold off;

ff = f; w = 2;
%ff([round(l/d)-w:round(l/d)+w,l-round(l/d)-w:l-round(l/d)+w])=0;
%ff([round(l/e)-w:round(l/e)+w,l-round(l/e)-w:l-round(l/e)+w])=0;

bb = abs(ifft(ff));
figure(4);
subplot(312);
plot(b); 
subplot(313);
plot(bb);ylim([0,1]);





% f= goertzel(b);%fft(b); %/length(a);   fwht
% figure(4);
% subplot(311);
% %plot(f);
% plot(real(f(1:round(length(a)/2)))); title(num2str(round(length(a)/e)));