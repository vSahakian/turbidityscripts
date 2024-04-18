%
%  Goal: Compare spectra from diferent years
%
clear all
close all
%
%  Define constants
%
day2min = 1440;
min2sec = 60;
day2sec = 1440*60;
tyears = [2016:1:2022];
nyears = length(tyears);
cc = linspecer(nyears);
iat = 1;
for jjj=2016:1:2022
    yearat = jjj;
    %
    %  Define where data are stored
    %  and load and define data
    %
    Tname = ['../Data_all_turbidity/reformated_BACAX_ntu_',num2str(yearat)];
    T_bark = readtable([Tname,'.csv']);
    timeA = T_bark.Var1;
    timeB = T_bark.Var2;
    data = T_bark.Var3;
    [y,m,d] = ymd(timeA);
    [h,mn,s] = hms(timeB);
    time_days = datenum(y,m,d,h,mn,s);
    time_min = day2min*time_days;
    time_sec = day2sec*time_days;
    delta_sec = time_sec(2) - time_sec(1)
    disp(['Now processing ',unique(y),' data.']);
    %
    %  Remove spikes using medfilt1
    %
    despike_len = 3;
    data_despike = medfilt1(data,despike_len);
    %
    %  Labels for data we will toss
    %
    %disp('Filtering!');
    %filt_low = 1e-6;
    %disp('&&&&& CAUTION! CHECK!!');
    filt_low = 0.0001
    %filt_low = 0.00001
    filt_high = 0.95*1/(2*60); % set to nyquest
    data_filt = eqfiltfilt(data_despike,filt_low,filt_high,delta_sec,4);
    %disp('Debug now');
    %keyboard
    % low freq signal data_filt = eqfiltfilt(data_despike,filt_low,filt_high,delta_sec,4);
    % Does not work data_filt = mtmspec(data_despike,delta_sec,4);
    % NW = 4;
    % nfft = 128;
    % data_filt =  mtmp_chatgbt(data_despike,NW,delta_sec,nfft,[filt_low filt_high]);
    lab_filt = [' filtbp: ',num2str(filt_low),'-',num2str(filt_high)];
    %
    %  Now convert time to mins relatiave to a given reference date
    %
    time_rel_days = time_days - min(time_days);
    time_rel_hours = time_rel_days*24;
    time_rel_min = time_rel_hours*60;

    seismogram = data;
    %time = time_days;
    time = time_sec;

    yyyy = unique(datestr(time_days,'yyyy'),'rows');
    disp(['Year ',yyyy,'; Read in ',num2str(length(time_days)),' data points']);
    disp(['      tmin = ',datestr(min(time_days))]);
    disp(['      tax = ',datestr(max(time_days))]);
    %disp('HACK!!!!***** Avoids gappy data!!!!');
    dt = unique(round(diff(time_rel_min(1:50))));
    disp(['      dt = ',num2str(dt),' in mins']);
    %
    % Compute nyquest
    %
    N = length(seismogram);
    Fs = 1 / (time(2) - time(1)); % Sampling frequency in seconds
    frequencies = Fs*(0:(N/2))/N; % Frequency axis

    seismogram_fft = fft(seismogram);
    amplitude_spectrum = 2*abs(seismogram_fft(1:N/2+1));

    %seismogram_fft_filt = fft(seismogram_filt);
    %amplitude_spectrum_filt = 2*abs(seismogram_fft_filt(1:N/2+1));

    % Plot the spectra
    hkeep(iat) = loglog(frequencies, amplitude_spectrum, 'b', 'LineWidth', 1,'Color',cc(iat,:));
    hold on
    %loglog(frequencies, amplitude_spectrum_filt, 'r', 'LineWidth', 1);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    %title(['Spectral Content of ',yyyy,' turbidity data, (Nsamples=',num2str(length(seismogram))]);
    title(['Spectral Content of turbidity data']);
    disp(['Length of data for year ',yyyy,' is ',num2str(length(time))]);
    grid on;

    hold on
    %yy = get(gca,'Ylim');
    %plot([filt_high filt_high],yy,'r--','LineWidth',2)
    %plot([filt_low filt_low],yy,'r--','LineWidth',2)
    %pname = [sub_dir,'plot_spectra_',num2str(yyyy),'.png'];
    %print(gcf,pname,'-dpng');
    iat = iat + 1;
end
set(gca,'YScale','log');
set(gca,'XScale','log');
set(gca,'FontSize',11);
legend(hkeep,['2016';'2017';'2018';'2019';'2020';'2021';'2022'],'Location','SouthWest');
pname = ['../plot_spectra_all_years.png']
print(gcf,pname,'-dpng');
