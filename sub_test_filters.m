function sub_test_filters(yearat,dirat)
%
%  Define constants
%
day2min = 1440;
min2sec = 60;
day2sec = 1440*60;
%
%  Define where data are stored
%
Tname = ['../Data_all_turbidity/reformated_BACAX_ntu_',num2str(yearat)];
%
%  Define subdirectories for plots and data
%
%Oname_data = [dirat,'data_',num2str(yearat),'/'];
%Oname_plots = [dirat,'plots_',num2str(yearat),'/'];
%system(['mkdir ',Oname_plots]);
%system(['mkdir ',Oname_data]);
%
%  Load data
%
T_bark = readtable([Tname,'.csv']);
timeA = T_bark.Var1;
timeB = T_bark.Var2;
data = T_bark.Var3;
[y,m,d] = ymd(timeA);
[h,mn,s] = hms(timeB);
time_days = datenum(y,m,d,h,mn,s);
time_min = day2min*time_days;
time_sec = day2sec*time_days;
delta_sec = time_sec(2) - time_sec(1);
%
%  Now convert time to mins relatiave to a given reference date
%
time_rel_days = time_days - min(time_days);
time_rel_hours = time_rel_days*24;
time_rel_min = time_rel_hours*60;
disp(['Read in ',num2str(length(time_days)),' data points']);
disp(['      tmin = ',datestr(min(time_days))]);
disp(['      tax = ',datestr(max(time_days))]);
disp('HACK!!!!***** Avoids gappy data!!!!');
dt = unique(round(diff(time_rel_min(1:50))));
disp(['Time Step = ',num2str(dt)]);
disp(['      dt = ',num2str(dt),' in mins']);
%disp(['Now processing ',unique(y),' data.']);
%
%  Remove spikes using medfilt1
%
despike_len = 3;
data_despike = medfilt1(data,despike_len);
labS = [' despiked using medfilt1 with Nlen=',num2str(despike_len)];
%
%  Filter data, testing different low frequency value
%
disp('Filtering!');
%filt_low = 1e-5; % Fav?
%filt_low = 1e-6;
%disp('&&&&& CAUTION! CHECK!!');
filt_low = 1e-3;
%filt_low = 0.0001
%filt_low = 0.00001
filt_high = 0.95*1/(2*60); % set to nyquest
filt_test = [1e-3, 1e-4, 1e-5, 1e-6];
filt_order = 4; 
%filt_order = 10; 
%filt_order = 50;  % Too high
%filt_order = 100; % Too high!
for ifilt = 1:length(filt_test)
    filt_low = filt_test(ifilt);
    lab_filt = [' filtbp: ',num2str(filt_low),'-',num2str(filt_high)];
    data_filt_old = eqfiltfilt(data_despike,filt_low,filt_high,delta_sec,4);
    D1 = designfilt('bandpassfir','FilterOrder',filt_order,...
        'CutoffFrequency1',filt_low,'CutoffFrequency2',filt_high,...
        'SampleRate',delta_sec);
    data_filt_new = filter(D1,data_despike);

    % 'arbmagfir'          Arbitrary magnitude FIR response
    % 'bandpassfir'        Bandpass FIR response
    % NO 'bandpassiir'        Bandpass IIR response
    % NO 'bandstopfir'        Bandstop FIR response
    % NO 'bandstopiir'        Bandstop IIR response
    % 'differentiatorfir'  Differentiator FIR response
    % NO 'highpassfir'        Highpass FIR response
    % NO 'highpassiir'        Highpass IIR response
    % 'hilbertfir'         Hilbert transformer FIR response
    % NO 'lowpassfir'         Lowpass FIR response
    % NO 'lowpassiir'         Lowpass IIR response
    %
    %  Title information
    %
    t2=[num2str(unique(y)),': filtbp: ',num2str(filt_low),'-',num2str(filt_high)];
    axis tight
    %
    % Plot raw data
    %
    figure
    plot(time_days,data,'k-','Color',[0.8 0.8 0.8]);
    hold on
    plot(time_days,data_filt_old,'b-');
    plot(time_days,data_filt_new,'g-');
    hold on
    axis tight
    title(['Data ',num2str(yearat),'; filtbp: ',num2str(filt_low),' to ',num2str(filt_high),' filtOrder=',num2str(filt_order)]);
    grid on
    ylabel('Turbidity (units)')
    xlabel('Time')
    set(gca,'FontSize',12)
    legend('Original Data','Filtered Data (old)','Filtered Data (new)');
    yy = get(gca,'Ylim');
    datetick('x','KeepLimits');
    set(gca,'Xlim',[datenum(2022,01,01,0,0,0) datenum(2022,01,05,0,0,0)]);
    %pname = [dirat,Oname_plots,'plot_data_',num2str(yearat),'.png'];
    %print(gcf,pname,'-dpng');
end
%
