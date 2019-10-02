%lag in ms
function [ pcor, timescale, train1smooth, train2smooth] = ...
                    timeCorrelationSmoothed( train1,train2,p)

train1smooth = train1; train2smooth = train2;
if isfield(p,'movmean') && p.movmean>0
    %'moving mean'
    %tic
    %win=hamming(p.hamwin);
    %train1smooth=conv(train1,win,'same');train2smooth=conv(train2,win,'same');
    %win = win/sum(win);
    train1smooth = movmean(train1,p.movmean);
    train2smooth = movmean(train2,p.movmean);
    %
    %toc
end
if isfield(p,'sigma') && p.sigma > 0
    'low pass(move mean)'
    tic
    train1smooth = lowPassFilter(train1smooth, lowf);
    train2smooth = lowPassFilter(train2smooth, lowf);
    toc
    %sgolayfilt(train1smooth,6,framelen)
end
if isfield(p,'bandf') && p.bandf >= 0
    'bandpass'
    %train1smooth = lowPassFilter(train1smooth, lowf);
    %train2smooth = lowPassFilter(train2smooth, lowf);
    train1smooth = bandPassFilter(train1smooth, p.bandf, p.bandw);
    train2smooth = bandPassFilter(train2smooth, p.bandf, p.bandw);
    %sgolayfilt(train1smooth,6,framelen)
end
%  figure, plot(train1Ham), hold on, plot(train1 + 1);
%correlation
%pcor = xcorr(train1Ham-mean(train1Ham), train2Ham-mean(train2Ham),lag,'coef');%SAME AS NEXT LINE
%p{i,3} = [pmx (pmxi - ((length(p{i,1})-1)/2) - 1)/1000]; % dont ask
%smoothing correlation
% if isfield(p,'lag') && p.sigma ~= 0
%     %[b,a] = butter(7,sigma); %LOW PASS 0.15%6th order, fc/fs/2 determined empiracally
%     %pcor = filtfilt(b,a,pcor);
%     %win=hamming(sigma);
%     %pcor=conv(pcor,win,'same');
%     
% end
%DO CORRELATION
if isfield(p,'lag') && p.lag > 0
    pcor = xcov(train1smooth, train2smooth,round(p.lag),'coef'); % TRY AS MATRIX 'coef'
else
    pcor = corr(train1smooth', train2smooth');
end
timescale=( (1:length(pcor)) - ((length(pcor)-1)/2) - 1)/1000;

end




%normalized frequency = freq *2 / sampling rate
%normalized frequency 0.1. Converted to Hz, this is 0.1*4000/2 = 200 Hz. 4000 sampling rate



function data = lowPassFilter(input, f)
%HIGH PASS FILTERING
% http://www.mathworks.com/help/dsp/ref/fdesign.bandpass.html
% All frequency values are in Hz.
% Construct an FDESIGN object and call its BUTTER method.
Fpass = f;      %600    % First Passband Frequency
Fstop = f+1;                % First Stopband Frequency
Apass  = .1;            % Passband Ripple (dB)
Astop = 100;            % First Stopband Attenuation (dB) %100
fs = 1000;     %32000  %data in ms
Fs = fs;
%toDecibal = 20*log(10);
%Hd = design(fdesign.lowpass(Fpass, Fstop, Apass, Astop, fs),'butter');
N = 8;
Hd = designfilt('lowpassiir', ...
    'FilterOrder',N, ...
    'PassbandFrequency',Fpass, ...
    'PassbandRipple',Apass,'StopbandAttenuation',Astop, ...
    'SampleRate',fs);
%fvtool(Hd);
%tic
data = filtfilt(Hd,input);
%toc
end

function data = bandPassFilter(input,bf,bw)
f1 = max(1e-9,bf-bw);
f2 = min(500-1e-9,bf+bw);
fs = 1000;%sampling frequency
data = bandpassfft(input,f1,f2,fs);

%HIGH PASS FILTERING
% http://www.mathworks.com/help/dsp/ref/fdesign.bandpass.html
% All frequency values are in Hz.
%{
f
Fpass1 = f-19;      %600    % First Passband Frequency
Fpass2 = f+0;                % First Stopband Frequency
Apass  = .1;            % Passband Ripple (dB)
Astop1 = 100;            % First Stopband Attenuation (dB) %100
Astop2 = 100;
Fs = 1000;     %32000  %data in ms
N = 2;
d = designfilt('bandpassiir', ...
  'FilterOrder',N, ...
  'PassbandFrequency1', Fpass1,'PassbandFrequency2', Fpass2, ...
  'StopbandAttenuation1', Astop1, 'PassbandRipple', Apass, ...
  'StopbandAttenuation2', Astop2, ...
  'SampleRate', Fs);

Hd = design(fdesign.bandpass(Fstop, Fpass, Astop, Apass, Fs),'butter');
data = filtfilt(d,input);
%}

end

function data = highPassFilter(input, f)
%HIGH PASS FILTERING
% http://www.mathworks.com/help/dsp/ref/fdesign.bandpass.html
% All frequency values are in Hz.
% Construct an FDESIGN object and call its BUTTER method.
f
Fstop =f-1;   %590    % First Stopband Frequency
Fpass =f;      %600    % First Passband Frequency This value is used because it interferes with neuron activity.
Astop = 10;            % First Stopband Attenuation (dB)
Apass  = 5;            % Passband Ripple (dB)
fs = 1000;     %32000  %data in ms
%toDecibal = 20*log(10);
Hd = design(fdesign.highpass(Fstop, Fpass, Astop, Apass, fs),'butter');
data = filtfilt(Hd.sosMatrix,Hd.ScaleValues,input);
end

