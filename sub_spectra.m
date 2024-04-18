function sub_spectra(time,seismogram,seismogram_filt,sub_dir)
% Generate a sample seismogram 
% Assuming 'seismogram' contains your seismogram data
% Assuming 'time' contains the corresponding time vector in days

% Compute the FFT

disp(['sub_dir=',sub_dir]);

%seismogram = trace;
%time = days_2_sec*sampletimes_days; % Assumes time is in days

day2sec = 1440*60
yyyy = unique(datestr(time,'yyyy'),'rows')
time = day2sec*time;

N = length(seismogram);
Fs = 1 / (time(2) - time(1)); % Sampling frequency
frequencies = Fs*(0:(N/2))/N; % Frequency axis

seismogram_fft = fft(seismogram);
amplitude_spectrum = 2*abs(seismogram_fft(1:N/2+1));

seismogram_fft_filt = fft(seismogram_filt);
amplitude_spectrum_filt = 2*abs(seismogram_fft_filt(1:N/2+1));

% Plot the spectra
figure;
loglog(frequencies, amplitude_spectrum, 'b', 'LineWidth', 2);
hold on
loglog(frequencies, amplitude_spectrum_filt, 'r', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title(['Spectral Content of ',yyyy,' turbidity data, (Nsamples=',num2str(length(seismogram))]);
grid on;

hold on
%yy = get(gca,'Ylim');
%plot([filt_high filt_high],yy,'r--','LineWidth',2)
%plot([filt_low filt_low],yy,'r--','LineWidth',2)
pname = [sub_dir,'plot_spectra_',num2str(yyyy),'.png'];
print(gcf,pname,'-dpng');
