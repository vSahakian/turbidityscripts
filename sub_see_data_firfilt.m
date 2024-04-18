function sub_see_data_firfilt(yearat,dirat)
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
%  Define directory where where results and metadata will be storred
%
%Mname_filt = [dirat,num2str(yearat),'_metadata_filtered.csv']
%fidM_filt = fopen(Mname_filt,'w+');

%Mname_nofilt = [dirat,num2str(yearat),'_metadata_unfiltered.csv']
%fidM_unfilt = fopen(Mname_nofilt,'w+');
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
%disp(['Now processing ',unique(y),' data.']);
%
%  Remove spikes using medfilt1
%
despike_len = 3;
data_despike = medfilt1(data,despike_len);
labS = [' despiked using medfilt1 with Nlen=',num2str(despike_len)];
%
%  Labels for data we will toss
%
%toss_lab = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'];
%
%  Plot oringal data and Filtered data
%
figure
clf
subplot(2,1,1)
plot(time_days,data,'b-')
hold on
plot(time_days,data_despike,'r-')
axis tight
datetick('x')
title({[num2str(yearat),' data: N = ',num2str(round(max(time_days)-min(time_days))),' days'],labS});
ylabel('Turbidity');
%
%  Filter data
%
disp('Filtering!');
filt_low = 1e-5; % Fav?
%filt_low = 1e-6;
%disp('&&&&& CAUTION! CHECK!!');
%filt_low = 1e-3;
%filt_low = 0.0001
%filt_low = 0.00001
filt_high = 0.95*1/(2*60); % set to nyquest
%data_filt = eqfiltfilt(data_despike,filt_low,filt_high,delta_sec,4);
data_filt = eqbandpassfir(data_despike,filt_low,filt_high,delta_sec,4);
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
disp(['Read in ',num2str(length(time_days)),' data points']);
disp(['      tmin = ',datestr(min(time_days))]);
disp(['      tax = ',datestr(max(time_days))]);
disp('HACK!!!!***** Avoids gappy data!!!!');
dt = unique(round(diff(time_rel_min(1:50))));
disp(['      dt = ',num2str(dt),' in mins']);
%
%  Flag to plot curve fit uncerts
%
ifitfig = 1;
%
%  Divide data into approximately weeks
%   & compute running average for each approximate month
%
tstep_week = 10081; % this is approximatly a week in minutes
tstep_week_days = tstep_week/(60*24); % convert from miinutes to days
%
%  Title information
%
t2=[num2str(unique(y)),': tstep=',num2str(round(tstep_week_days)),'days, filtbp: ',num2str(filt_low),'-',num2str(filt_high)];
axis tight
%
%  Compute moving mean in order to establish the noise floor
%
%MMweek = movmean(data_filt,tstep_week); % these are arrays
%SSweek = movstd(data_filt,tstep_week);
%TThresh_week = MMweek + 2*SSweek;
%
% Plot raw data
%
subplot(2,1,2)
plot(time_days,data,'k-','Color',[0.8 0.8 0.8]);
hold on
plot(time_days,data_filt,'g-','Color',[0.2 0.8 0.4]);
hold on
axis tight
title(['Data ',num2str(yearat),'; filtbp: ',num2str(filt_low),'-',num2str(filt_high),' original data==grey, filtered==green'])
grid on
ylabel('Turbidity (units)')
xlabel('Time')
set(gca,'FontSize',12)
yy = get(gca,'Ylim');
%plot(time_days,MMweek,'k-','Color',[0.8 0.6 0.4],'LineWidth',2);
datetick('x');
axis tight
%legend('Original Data','Filtered Data',['Running mean of consecutive points (N=',num2str(round(tstep_week_days)),' days)'],'Location','southoutside');
%legend('Original Data','Filtered Data','Location','southoutside');
set(gca,'FontSize',11);
%pname = [dirat,Oname_plots,'plot_data_',num2str(yearat),'.png'];
%print(gcf,pname,'-dpng');
%
