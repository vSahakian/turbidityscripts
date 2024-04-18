%
%  Goal: Plot all years of data on a single plot
%
%  Define constants
%
day2min = 1440;
min2sec = 60;
day2sec = 1440*60;
%
%  Load all data 
%
years = [2016; 2017; 2018; 2019; 2020; 2021; 2022];
nyears = length(years);

time_all = [];
data_all = [];
for jloop = 1:nyears
    yearat = years(jloop);
    Tname = ['../Data_all_turbidity/reformated_BACAX_ntu_',num2str(yearat)];
    T_bark = readtable([Tname,'.csv']);

    timeA = T_bark.Var1;
    timeB = T_bark.Var2;
    data = T_bark.Var3;
    disp([num2str(yearat),': N=',num2str(length(data))]);

    [y,m,d] = ymd(timeA);
    [h,mn,s] = hms(timeB);
    time_days = datenum(y,m,d,h,mn,s);
    npoints(jloop) = length(time_days);
    %
    time_all = [time_all; time_days];
    data_all = [data_all; data];
    eval(['time_',num2str(yearat),'=time_days;']);
    eval(['data_',num2str(yearat),'=data;']);

end % End jloop

figure(1)
clf
plot(time_all,data_all,'g.-');
hold on
axis tight
datetick('x');
yy = get(gca,'Ylim');
for i=1:nyears
   plot([datenum(years(i),01,01,0,0,0) datenum(years(i),01,01,0,0,0)],yy,'b--');
   h=text(datenum(years(i),01,10,0,0,0),yy(2)-0.03*(yy(2)-yy(1)),[num2str(years(i)),': ',num2str(round(npoints(i)./1440)),' days']);
   set(h,'FontSize',7);
end
print(gcf,'plot_data_all.png','-dpng');

figure(2)
clf
orient tall
subplot(nyears,1,1);
plot(time_2016,data_2016);
set(gca,'Xlim',[datenum(2016,01,01,0,0,0) datenum(2016,12,31,59,59,59)]);
datetick('x','Keeplimits');
title(['2016, N=',num2str(length(time_2016))]);
ylabel('Turbidity (NTU)');
%set(gca,'FontSize',12);
subplot(nyears,1,2);
plot(time_2017,data_2017);
title(['2017, N=',num2str(length(time_2017))]);
set(gca,'Xlim',[datenum(2017,01,01,0,0,0) datenum(2017,12,31,59,59,59)]);
datetick('x','Keeplimits');
ylabel('Turbidity (NTU)');
%set(gca,'FontSize',12);
subplot(nyears,1,3);
plot(time_2018,data_2018);
title(['2018, N=',num2str(length(time_2018))]);
set(gca,'Xlim',[datenum(2018,01,01,0,0,0) datenum(2018,12,31,59,59,59)]);
datetick('x','Keeplimits');
ylabel('Turbidity (NTU)');
%set(gca,'FontSize',12);
subplot(nyears,1,4);
plot(time_2019,data_2019);
title(['2019, N=',num2str(length(time_2019))]);
set(gca,'Xlim',[datenum(2019,01,01,0,0,0) datenum(2019,12,31,59,59,59)]);
datetick('x','Keeplimits');
ylabel('Turbidity (NTU)');
%set(gca,'FontSize',12);
subplot(nyears,1,5);
plot(time_2020,data_2020);
title(['2020, N=',num2str(length(time_2020))]);
set(gca,'Xlim',[datenum(2020,01,01,0,0,0) datenum(2020,12,31,59,59,59)]);
datetick('x','Keeplimits');
ylabel('Turbidity (NTU)');
%set(gca,'FontSize',12);
subplot(nyears,1,6);
plot(time_2021,data_2021);
title(['2021, N=',num2str(length(time_2021))]);
datetick('x','Keeplimits');
set(gca,'Xlim',[datenum(2021,01,01,0,0,0) datenum(2021,12,31,59,59,59)]);
ylabel('Turbidity (NTU)');
%set(gca,'FontSize',12);
subplot(nyears,1,7);
plot(time_2022,data_2022);
set(gca,'Xlim',[datenum(2022,01,01,0,0,0) datenum(2022,12,31,59,59,59)]);
datetick('x','Keeplimits');
title(['2022, N=',num2str(length(time_2022))]);
ylabel('Turbidity (NTU)');
%set(gca,'FontSize',12);
print(gcf,'plot_data_by_year.png','-dpng');



