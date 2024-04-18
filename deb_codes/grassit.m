%
%  Goal: Determine outlier thresholds for the turbidity data
%
clear all
%
%  Load 3 different datasets
%    1. M7+ global
%    2. M5+ regional (within 1000 km)
%    3. M3+ local (within 300 km)
%
%load ../Data_earthquakes/eqs_over3_xyzm.dat
%load ../Data_earthquakes/eqs_over5_xyzm.dat
%load ../Data_earthquakes/eqs_over7_xyzm.dat
%
load ../Data_earthquakes/earthquakes_and_strain.mat

disp('Combinine M3+, M5+ and M7+ data into a single array');

junk3 = s3.lat3;
junk5 = s5.lat5;
junk7 = s7.lat7;

n3 = length(junk3)
n5 = length(junk5)
n7 = length(junk7)

eq_lat = extractfield(s3,'lat3');
eq_lon = extractfield(s3,'lon3');
eq_dep = extractfield(s3,'dep3');
eq_mag = extractfield(s3,'mag3');
eq_dist = extractfield(s3,'dist3');
eq_time = extractfield(s3,'time3');
eq_strain = extractfield(s3,'strain3');

eq_lat = [eq_lat, extractfield(s5,'lat5')];
eq_lon = [eq_lon, extractfield(s5,'lon5')];
eq_dep = [eq_dep, extractfield(s5,'dep5')];
eq_mag = [eq_mag, extractfield(s5,'mag5')];
eq_dist = [eq_dist, extractfield(s5,'dist5')];
eq_time = [eq_time, extractfield(s5,'time5')];
eq_strain = [eq_strain, extractfield(s5,'strain5')];

eq_lat = [eq_lat, extractfield(s7,'lat7')];
eq_lon = [eq_lon, extractfield(s7,'lon7')];
eq_dep = [eq_dep, extractfield(s7,'dep7')];
eq_mag = [eq_mag, extractfield(s7,'mag7')];
eq_dist = [eq_dist, extractfield(s7,'dist7')];
eq_time = [eq_time, extractfield(s7,'time7')];
eq_strain = [eq_strain, extractfield(s7,'strain7')];
%
%  Define time window of earthquake catalog
%
t1_year = datestr(min(eq_time),'yyyy');
t1_mo = datestr(min(eq_time),'mm');

t2_year = datestr(max(eq_time),'yyyy');
t2_mo = datestr(max(eq_time),'mm');
tlab = [t1_year,'-',t1_mo,'-01 - ',t2_year,'-',t2_mo,'-30'];
%
%  Define constants
%
day2min = 1440;
day2sec = 1440*60
%
%  Count number of 3+ local, 5+ regional and 7+ global
%
i3 = find(eq_mag<5);
i5 = intersect(find(eq_mag>=5),find(eq_mag<7));
i7 = find(eq_mag>=7);
%
%  Selected earthquakes with high strains
%
strain_thresh = 0.05;
ipick = find(eq_strain>=strain_thresh)
%[junk,iii]  = sort(eq_strain);
%ipick = iii(1:10);
for i=1:length(ipick)
    disp([datestr(eq_time(ipick(i))),' M',num2str(eq_mag(ipick(i))),' ',num2str(eq_strain(ipick(i)))]);
end
nfav = length(ipick);
disp(['Number of high strain events: ',num2str(nfav)]);
%
% Strain plot
%
figure(3)
clf
semilogy(eq_mag(i7),eq_strain(i7),'go','Color',[0.0 0.4 0.8],'MarkerFaceColor',[0.2 0.6 0.8],'MarkerSize',7);
hold on
semilogy(eq_mag(i5),eq_strain(i5),'ro','Color',[0.8 0.6 1.0],'MarkerFaceColor',[0.8 0.4 1.0],'MarkerSize',7);
semilogy(eq_mag(i3),eq_strain(i3),'yo','Color',[0.8 0.7 0.5],'MarkerFaceColor',[0.9 0.8 0.6],'MarkerSize',7);
xx = get(gca,'Xlim');
set(gca,'Xlim',[3.0 xx(2)]);
yy = get(gca,'Ylim');
%set(gca,'Ylim',[yy(1) 1.6]);
%yy = get(gca,'Ylim');
plot(xx,[strain_thresh strain_thresh],'k--');
title(['Assumed dynamic strain threshold ',num2str(strain_thresh),' (dashed line)']);
%plot(xx,[1 1],'k:');
xlabel('Earthquake Magnitude');
ylabel('Strain');
set(gca,'FontName','Helvetica','FontSize',12)
lab1 = ['Global earthquakes (M7+, ',tlab,', N=',num2str(n7),')'];
lab2 = ['Regional earthquakes within 1000 km (M5-7, ',tlab,', N=',num2str(n5),')'];
lab3 =  ['Local earthquakes within 300 km (M3-5, ',tlab,', N=',num2str(n3),')'];
%h=legend(['Regional earthquakes within 300 km (M5+, N=',num2str(n5),')'],['Global earthquakes (M7+, N=',num2str(n7),')'],'Location','SouthEast');
h=legend(lab1,lab2,lab3,'Location','SouthEast');
set(h,'FontSize',8)
set(gca,'FontName','Helvetica','FontSize',12);
text(xx(1)+0.01*(xx(2)-xx(1)),1.0,'(a)');
set(gca,'Yscale','log')
grid on
%
%  Define plot directory
%
pdir = 'plots';
if ~exist(pdir,'dir')
    system('mkdir plots');
end
%
%  Load all avaiable turbidity data
%
years_want = [2016:1:2022];
turbid_orig = [];
time_turbid_days_orig = [];
for jj=1:length(years_want)
    yearat = years_want(jj);
    clear Tname y m d h mn s
    Tname = ['../Data_all_turbidity/reformated_BACAX_ntu_',num2str(yearat)];
    T_bark = readtable([Tname,'.csv']);
    timeA = T_bark.Var1;
    timeB = T_bark.Var2;
    [y,m,d] = ymd(timeA);
    [h,mn,s] = hms(timeB);
    time_turbid_days_orig = [time_turbid_days_orig; datenum(y,m,d,h,mn,s)];
    turbid_orig = [turbid_orig; T_bark.Var3];
end % end jj loop over years of interest
disp('Loaded turbidity data');
%
figure(10)
clf
plot(time_turbid_days_orig,turbid_orig,'b.');
hold on
plot(eq_time(ipick),zeros(size(eq_time(ipick))),'rp','MarkerSize',18,'MarkerFaceColor',[0.8 0.2 0.1]);
npick = length(ipick);
for i=1:npick
   text(eq_time(ipick(i)),0,num2str(eq_strain(ipick(i)),2));
end
title('Turbidity data at BACX, high strain events (starts) and strains listed');
datetick('x','yyyy');
print(gcf,'-dpng',[pdir,'/plot_turb_data_avail.png']);
%
%  Currently, only have turbdidty data up to '03-Feb-2020 10:47:30'
%
itoss=[];
for kk=1:length(ipick)
    gaptest(kk) = min(abs(eq_time(ipick(kk))-time_turbid_days_orig))
    if gaptest(kk)>1
       disp(['Need to toss high strain event: ',datestr(eq_time(kk))]);
       itoss = [itoss;kk];
    end
end
%years_pick = str2num(datestr(eq_time(ipick),'yyyy'));
%itoss = find(years_pick>=2020);
disp('****** CAUTION!  *****');
disp('****** CAUTION!  *****');
disp('****** CAUTION!  *****');
disp('****** CAUTION!  *****');
disp('****** CAUTION!  *****');
disp(['Need to toss ',num2str(length(itoss)),' events that do not have associated turbidity data']);
ipick(itoss) = [];
%
%  Deinfe static parameters
%
%Pwant = [ 95 99];
%Pwant = [ 95 98 97];
Pwant = [ 95 ];
%Pwant = [ 90 ];
%Pwant = [ 40 ];
%ttest_hours = [  12 8 ]; % Caution, negative numbers if using 2-hours ... sometimes
%ttest_hours = [  12 8 ]; % Caution, negative numbers if using 2-hours ... sometimes
ttest_hours = [ 12 ];
tmoveit_hours = 1;
ttest_days = ttest_hours./24;
tmovit_days = tmoveit_hours./24;
%
%  Load turbidity data
%
npick = length(ipick);
disp(['About to process ',num2str(npick),' high strain events']);
for iii=1:npick
%for iii=4:4
    yearat = datestr(eq_time(ipick(iii)),'yyyy')
    time_main = eq_time(ipick(iii));

    %Tname = ['../Data_all_turbidity/reformated_BACAX_ntu_',num2str(yearat)];
    %T_bark = readtable([Tname,'.csv']);
    %timeA = T_bark.Var1;
    %timeB = T_bark.Var2;
    %turbid_orig = T_bark.Var3;
    %[y,m,d] = ymd(timeA);
    %[h,mn,s] = hms(timeB);
    %time_turbid_days_orig = datenum(y,m,d,h,mn,s);
    %disp(['Now processing ',datestr(eq_time(ipick(iii))),' data.']);
    %
    %   Define time window
    %
%nmonths = 2;
nmonths = 1;
    disp(['Limiting turbidity time window to plus/minus ',num2str(nmonths),' months']);
    timeA = eq_time(ipick(iii)) - nmonths*30.4;
    timeB = eq_time(ipick(iii)) + nmonths*30.4;
    ikeepA = find(time_turbid_days_orig>=timeA);
    ikeepB = find(time_turbid_days_orig<=timeB);
    ikeep = intersect(ikeepA,ikeepB);
    time_turbid_days = time_turbid_days_orig(ikeep);
    turbid = turbid_orig(ikeep);
    %time_turbid_min = day2min*time_turbid_days;
    %time_turbid_sec = day2sec*time_turbid_days;
    delta_turbid_sec = s(2) - s(1);

    ieq1 = find(eq_time>=timeA);
    ieq2 = find(eq_time<=timeB);
    ieq = intersect(ieq1,ieq2);
    time_eqs = eq_time(ieq);
disp('Update .... confusing!');
    neqs = length(time_eqs);
    disp(['Number of earthquakes within this 4 month time period: ',num2str(neqs)]);
    %
    %  Loop over all percentages
    %
    for nnn = 1:length(Pwant);
        %
        %  Define what percentage to use in the percentile tests
        %
        Pat  = Pwant(nnn);
        disp(['Working on ',num2str(Pat),'%']);
        %
        %  Clear all tsave arrays
        %
      tsave_outliers_24 = [];
      tsave_outliers_12 = [];
      tsave_outliers_8 = [];
      tsave_outliers_5 = [];
      tsave_outliers_2 = [];
      tsave_outliers_nfound_24 = [];
      tsave_outliers_nfound_12 = [];
      tsave_outliers_nfound_8 = [];
      tsave_outliers_nfound_5 = [];
      tsave_outliers_nfound_2 = [];
      tsave_thresh_24 = [];
      tsave_thresh_12 = [];
      tsave_thresh_8 = [];
      tsave_thresh_5 = [];
      tsave_thresh_2 = [];
      tsave_percentile_24 = [];
      tsave_percentile_12 = [];
      tsave_percentile_8 = [];
      tsave_percentile_5 = [];
      tsave_percentile_2 = [];
      tsave_yr_min = [];
      tsave_yr_max = [];
        %
        %  define time and turbid data
        %
        Dname = 'Fabio';
        Pname_base = 'plot_Fabio_proposed';
        Sname = 'Fabio_outlie';

        disp(['DEBUG: Number of time_turbid_days points: ',num2str(length(time_turbid_days))]);
        time_temp = time_turbid_days;
        turbid_temp = turbid;
        %
        %  Define stuff
        %
        N = length(time_turbid_days);
        %disp('Using all data, no need to loop over 3-year time windows ....');
        %
        % Make sure data is time sorted
        %
        N = length(time_temp);
        disp(['DEBUG: Number of time_temp points: ',num2str(N)]);
        [time_sort,isort] = sort(time_temp);
        turbid_sort = turbid_temp(isort);
        N = length(time_sort);
        time = time_sort;
        turid = turbid_sort;
        N = length(time);
        disp(['DEBUG: Number of points: ',num2str(N)]);
        %
        %  Define time window of interest
        %
        tstart = min(time);
        tend = max(time);
        %
        %  Loop over all test window durations and assign a moving window duration
        %
        for i=1:length(ttest_days)
            jat = i;
            testwin_days =ttest_days(i);
            testwin_hours =testwin_days*24;
            disp(['testwin_hours=',num2str(testwin_hours)]);
            t1 = min(time);
            t2 = t1+testwin_days;
            iat = 1;
            clear nfound timestamp
disp('DEBUG: Start while loop');
            while(t2<tend)
                iselect = find(and(time<t2,time>=t1)); % identify data snippet of interest
                nfound(iat) = sum(turbid(iselect));
                timestamp(iat) = t1;
                %timestamp(iat) = t2;
                %timestamp(iat) = t1+0.5*(t2-t1);
                t1 = t1+tmovit_days;
                t2 = t1+testwin_days;
                iat = iat + 1;
            end
    %keyboard
            %disp(['DEBUG iat=',num2str(iat),' ttest_hours=',num2str(ttest_hours(i))]);
            %disp(['DEBUG tmovit_days=',num2str(tmovit_days),' testwin_days=',num2str(testwin_days)]);

            %ioutlie = sub_cheb(nfound,Nstd);
            %[ioutlie,nthresh] = sub_percentile(nfound,Pat);
            if testwin_hours<3
                disp('CAUTION: Tail extreme, using log!!!');
                [ioutlie,nthresh] = sub_percentile_turbid(log(nfound),Pat);
            else
                [ioutlie,nthresh] = sub_percentile_turbid(nfound,Pat);
            end
            noutlier = length(ioutlie);
            %
            %if jat==1
            %   subplot(4,1,1)
            %   hist(nfound,[min(nfound):1:max(nfound)])
            %   xlabel('Nfound')
            %   ylabel('Number')
            %   %title(['Quick debug histogram of twin=24 nfound, NSTD=',num2str(Nstd)]);
            %   hold on
            %   yy = get(gca,'Ylim');
            %   %plot([stat_mean stat_mean],yy,'r--')
            %   %h=plot([stat_thresh stat_thresh],yy,'r-')
            %   %set(h,'LineWidth',2);
            %   %legend('Hist','Mean','Thresh','Location','NorthEastOutside');
            %   plot([nthresh nthresh],yy,'r--');
            %   legend('Hist','Percentile Threshold','Location','NorthEastOutside');
            %end

            figure(2)
            clf
            subplot(2,1,1)
            if testwin_hours<3
                h=plot(timestamp,log(nfound),'k.','Color',[0.0 0.6 0.4]);
                hold on
                axis tight
                xx = get(gca,'Xlim');
                plot(timestamp(ioutlie),log(nfound(ioutlie)),'ro','Color',[1.0 0.6 0.0],'LineWidth',1)
                ylabel('log(Turbidity)')
                plot(xx,[log(nthresh) log(nthresh)],'--','Color',[1.0 0.6 0.0],'LineWidth',2);
            else
                %timeA = datenum(2018,10,29,18,58,36);
                %timeB = datenum(2018,11,15,10,32,24);
                %pp=patch([timeA timeB timeB timeA timeA],[0 0 83.0 83.0 0],'y');
                %set(pp,'FaceAlpha',0.50,'LineStyle','none');
                hold on
                h=plot(timestamp,(nfound),'k.','Color',[0.0 0.6 0.4]);
                hold on
                axis tight
                xx = get(gca,'Xlim');
                plot(timestamp(ioutlie),(nfound(ioutlie)),'ro','Color',[1.0 0.6 0.0],'LineWidth',1)
                ylabel('Turbidity')
                %ylabel('S_T')
                plot(xx,[nthresh nthresh],'--','Color',[1.0 0.6 0.0],'LineWidth',2);
            end
            xlabel('Time')
            yy = get(gca,'Ylim');
            for i=1:neqs
                plot([time_eqs(i) time_eqs(i)],yy,'k:','Color',[0.8 0.8 0.8],'LineWidth',2);
            end
            plot([time_main time_main],yy,'r-','Color',[0.6 0.6 0.4],'LineWidth',3);

            %h=text( 737474.7311,      69.53555556, 'Barkley Canyon Axis');
            %set(h,'FontName','Helvetica','FontSize',14);

            %set(gca,'Xlim',[tstart tend]);
            datetick('x',28,'keeplimits')

            tlab2 =  ['Flagged ',num2str(length(ioutlie)), ...
                ' of ',num2str(iat),' windows (',num2str(100*length(ioutlie)./iat,2),'%)'];
            %tlab3 = ['Nthresh = ',num2str(nthresh),' Percent=',num2str(Pat),'%'];
            %tlab3 = ['N=',num2str(iat),', S_{thresh}=',num2str(nthresh),', ',num2str(Pat),'%'];
            tlab3 = ['N=',num2str(iat),', S_{thresh}=',num2str(nthresh)];

            %title({[Dname,': ',num2str(testwin_hours),' hr window; ',num2str(tmoveit_hours),' hr moving shift; Stats=Percentile ',tlab2],tlab3});
            %title({['\Deltat=',num2str(testwin_hours),' hrs; ',num2str(tmoveit_hours),' hr moving shift; ',tlab2],tlab3});
            %title(['\Deltat=',num2str(testwin_hours),' hrs; ',num2str(tmoveit_hours),' hr moving shift; ',tlab3]);
            %title({['\Deltatime=',num2str(testwin_hours),'h, moving shift=',num2str(tmoveit_hours),'h, ',tlab3],loc});
            %title(['Barkley Canyon Axis (',num2str(Pat),'%): \Deltatime=',num2str(testwin_hours),'h, moving shift=',num2str(tmoveit_hours),'h, ',tlab3]);
            title([datestr(time_main),' M',num2str(eq_mag(ipick(iii))),' strain=',num2str(eq_strain(ipick(iii)),2),' (',num2str(Pat),'%): \Deltatime=',num2str(testwin_hours),'h, moving shift=',num2str(tmoveit_hours),'h, ',tlab3]);
            %tout = [pdir,'/',Pname_base,'_paper.png'];
            %tout = [pdir,'/',Pname_base,'_Dt_',num2str(testwin_hours),'_',datestr(min(time),30),'_',datestr(max(time),30),'_Percent_',num2str(Pat),'.png'];
            tout = [pdir,'/plot_',num2str(iii),'.png'];
            set(gca,'FontSize',10)
            print(gcf,'-dpng',tout);
            %
            %   Save information
            %
            if (testwin_hours==24)
                tsave_outliers = timestamp(ioutlie);
                tsave_outliers_24 = [tsave_outliers_24, tsave_outliers];
                tsave_outliers_nfound_24 = [tsave_outliers_nfound_24, nfound(ioutlie)];
                tsave_thresh_24 = [tsave_thresh_24, nthresh];
                tsave_percentile_24 = [tsave_percentile_24, Pat];
            end
            if (testwin_hours==12)
                tsave_outliers = timestamp(ioutlie);
                tsave_outliers_12 = [tsave_outliers_12, tsave_outliers];
                tsave_outliers_nfound_12 = [tsave_outliers_nfound_12, nfound(ioutlie)];
                tsave_percentile_12 = [tsave_percentile_12, Pat];
                tsave_thresh_12 = [tsave_thresh_12, nthresh];
                %tsave_mean_12 = [tsave_mean_12, stat_mean];
                %tsave_std_12 = [tsave_std_12, stat_std];
            end
            if (testwin_hours==8)
                tsave_outliers = timestamp(ioutlie);
                tsave_outliers_8 = [tsave_outliers_8, tsave_outliers];
                tsave_outliers_nfound_8 = [tsave_outliers_nfound_8, nfound(ioutlie)];
                tsave_thresh_8 = [tsave_thresh_8, nthresh];
                tsave_percentile_8 = [tsave_percentile_8, Pat];
            end
            if (testwin_hours==5)
                tsave_outliers = timestamp(ioutlie);
                tsave_outliers_5 = [tsave_outliers_5, tsave_outliers];
                tsave_outliers_nfound_5 = [tsave_outliers_nfound_5, nfound(ioutlie)];
                tsave_thresh_5 = [tsave_thresh_5, nthresh];
                tsave_percentile_5 = [tsave_percentile_5, Pat];
            end
            if (testwin_hours==2)
                tsave_outliers = timestamp(ioutlie);
                tsave_outliers_2 = [tsave_outliers_2, tsave_outliers];
                tsave_outliers_nfound_2 = [tsave_outliers_nfound_2, nfound(ioutlie)];
                tsave_thresh_2 = [tsave_thresh_2, nthresh];
                tsave_percentile_2 = [tsave_percentile_2, Pat];
            end
        end % end loop over the time-windows
        %
        %  Now save the outlier information for use in the paper_matchit code
        %
        if ismember(ttest_hours,24)
            save([Sname,'_percent',num2str(Pat),'_24.mat'],'tsave_outliers_24','tsave_outliers_nfound_24','tsave_thresh_24',...
            	'tsave_yr_min','tsave_yr_max','Dname');
        end
        if ismember(ttest_hours,12)
            save([Sname,'_percent',num2str(Pat),'_12.mat'],'tsave_outliers_12','tsave_outliers_nfound_12','tsave_thresh_12',...
                'tsave_yr_min','tsave_yr_max','Dname');
        end
        if ismember(ttest_hours,5)
            save([Sname,'_percent',num2str(Pat),'_5.mat'],'tsave_outliers_5','tsave_outliers_nfound_5','tsave_thresh_5',...
                'tsave_yr_min','tsave_yr_max','Dname');
        end
        if ismember(ttest_hours,8)
            save([Sname,'_percent',num2str(Pat),'_8.mat'],'tsave_outliers_8','tsave_outliers_nfound_8','tsave_thresh_8',...
                'tsave_yr_min','tsave_yr_max','Dname');
        end
        if ismember(ttest_hours,2)
            save([Sname,'_percent',num2str(Pat),'_2.mat'],'tsave_outliers_2','tsave_outliers_nfound_2','tsave_thresh_2',...
                'tsave_yr_min','tsave_yr_max','Dname');
        end
    end % end loop percentages
end % end loop over high strain events
