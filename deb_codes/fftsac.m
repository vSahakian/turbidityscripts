function [spec, phase, f] = fftsac(data, delta) ;

%Perform an FFT on a random signal.

y = fft(data,length(data));
phase=unwrap(angle(y));
spec = abs(y);          % or sqrt(y.*conj(y));
f = ((1/(length(data)*delta))*(0:length(y)/2))';
spec = delta*spec(1:length(f));
phase = phase(1:length(f));

% plot(f, phase*180/pi)
% title('Phase, unwrapped')
% xlabel('frequency')
% pause
% semilogy(f,spec)
% title('Amplitude')
% xlabel('frequency')
% pause
% keyboard
