%
%  Goal: Identify data gaps
%
%  Define constants
%
day2min = 1440;
min2sec = 60;
day2sec = 1440*60;

min2day = 0.000694444;
%
%  Loop over all years of data
%
years = [2016; 2017; 2018; 2019; 2020; 2021; 2022];
nyears = length(years);

cc_colors = brewermap(7,'Accent')
figure(10)
clf
for jloop = 1:nyears
    yearat = years(jloop);
    %T_bark = readtable('bark.csv');
    %Tname = 'reformated_BACAX_ntu_2019'; % Save name for later use
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
    delta_sec = time_sec(2) - time_sec(1);
    %
    % Indentify gaps
    %
    igap = find(diff(time_min)>1.1);
    ngap = length(igap);
    disp(['Year ',num2str(yearat),' has ',num2str(ngap),' gaps'])
    for i=1:ngap
        time_gap(i) = time_days(igap(i));
    end
    %
    % Write gap information to a file
    %
    fid = fopen(['data_gaps_',num2str(yearat),'.dat'],'w');
    dt_min_tot = 0;
    for i=1:ngap
        figure(9)
        clf
        plot(time_days,data,'kx-');
        hold on
        t1 = time_gap(i)-0.1;
        t2 = time_gap(i)+0.1;
        dt_day = time_days(igap(i)+1) - time_days(igap(i));
        dt_min = dt_day*1440;
dt_min_tot = dt_min_tot + dt_min;
        set(gca,'Xlim',[t1 t2]);
        yy = get(gca,'Ylim');
        plot([time_gap(i) time_gap(i)],yy,'r-');
        grid on
        datetick('x','KeepLimits');
        title(['Gap #',num2str(i),' of 7: ',datestr(time_gap(i),30),' gap=',num2str(dt_min),' mins']);
        set(gca,'FontSize',12);
        pname = ['plot_',num2str(yearat),'_gap_',num2str(i),'.png'];
        print(gcf,pname,'-dpng');
        %disp('Hit any key to continue');
        %pause
        %
        %  Write out the gap information
        %
        i1 = igap(i);
        i2 = i1+1;
        fprintf(fid,'%s\n',[datestr(time_days(i1),30),' ',datestr(time_days(i2),30),' ',num2str(time_min(i2)-time_min(i1))]);
    end
    fclose(fid);
    %
    % check start and end
    %
    if datenum(yearat,01,01,0,0,30)==time_days(1)
       imatch_start =  1;
    else
       disp('Start Time not in January');
    end
    if datenum(yearat,12,31,23,59,30)==time_days(end)
       imatch_end =  1;
    else
       disp('End Time not in December');
    end
disp([num2str(yearat),': has ',num2str(dt_min_tot),' mins of missing data (',datestr(time_days(1)),' - ',datestr(time_days(end))]);
    figure(10)
    subplot(2,1,1)
    plot(time_days,data,'ko-','Color',cc_colors(jloop,:),'MarkerSize',4)
    hold on

    %
    % Compute number of days of data
    %
    trange = time_days(end) - time_days(1) - min2day*dt_min_tot;
    if leapyear(yearat)
       timeper = 100*(trange/366);
    else
       timeper = 100*(trange/365);
    end
    disp([num2str(yearat),': has ',num2str(round(trange)),' days of data. ~',num2str(round(timeper)),'% of the year']);

end % End jloop
figure(10)
datetick('x','KeepLimits');
xlabel('Date');
ylabel('Turbidity');
grid on
title('ONC turbidity data (color by year)');
print plot_time_all_years.png -dpng
