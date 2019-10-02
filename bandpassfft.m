function o = bandpassfft(x, f1, f2, fs)
%remove mean
y = fft(x);%fft(x-mean(x));             % MEAN
%t = 1/fs;               % sampling period
l = length(x);       % length of signal
%t = (0:l-1)*t;        % time vector
%f = fs*(0:(l/2))/l;  % actual frequency, not normalized
% plot(f(1:1000),real(y(l:l+1000)));
%add mean

f1i = f1*l/fs;
f2i = f2*l/fs;
i = ones(1,l);
i(f1i:f2i) = 0; i(end-f2i:end-f1i) = 0; i = boolean(i);
iy = y; iy(i)=0;
o = ifft(iy);%+mean(x); %MEAN
figure(99);
a = subplot(211);
plotfft(a,x,fs);
a = subplot(212);
plotfft(a,o,fs);




end