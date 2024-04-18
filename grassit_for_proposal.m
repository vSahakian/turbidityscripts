%
%  Goal: Determine outlier thresholds for the turbidity data
%
clear all
%
%  Select dataset
%
iwant_2014 = 1;
iwant_2018_all = 0;
iwant_2018_some = 0;
iwant_2012 = 0;
%
%  Define plot directory
%
pdir = 'plots';
if ~exist(pdir,'dir')
   system('mkdir plots');
end
%
%
%  Deinfe static parameters
%Pwant = [ 95 99];
%Pwant = [ 95 98 97];
Pwant = [ 95 ];
%Pwant = [ 90 ];
%Pwant = [ 40 ];
%ttest_hours = [  12 8 ]; % Caution, negative numbers if using 2-hours ... sometimes
ttest_hours = [  12 8 ]; % Caution, negative numbers if using 2-hours ... sometimes
%ttest_hours = [  5 ];
%disp('CAUTION only 5 hours');
%ttest_hours = [  24 ];
%disp('CAUTION only 24 hours');
%ttest_hours = [  2 ];
%disp('CAUTION only 2 hours');
%ttest_hours = [  12 ];
%disp('CAUTION only 12 hours');
%ttest_hours = [  8 ];
%disp('CAUTION only 8 hours');
%ttest_hours = [  48 ];
%disp('CAUTION only 48 hours');
tmoveit_hours = 1;
ttest_days = ttest_hours./24;
tmovit_days = tmoveit_hours./24;
%
%  Load data
%
if iwant_2012==1
   load ../ONC_data/20120101T000523Z_20140101T000443Z.mat
   time_orig = data.time;
   turbid_orig = data.dat;
   time = time_orig;
   turbid = turbid_orig;
   iwantA = find(time_orig<datenum(2013,01,15,0,0,0));
   iwantB = find(time_orig>datenum(2012,09,01,0,0,0));
   iwant = intersect(iwantA,iwantB);
   time = time_orig(iwant);
   turbid = turbid_orig(iwant);
   met_loc = metadata.location;
   loc = [met_loc.locationName,' ',met_loc.stationCode];
end
if iwant_2014==1
   T_fabio = readtable('fabio_2014.csv');
   loc = ['Fabio: Barkley Canyon MidEast (May2014-Jan2015)'];
   timeX = T_fabio.Time;
   temp_cell = T_fabio.Temperature;
   sigma = T_fabio.Sigma_t_kg_m3_;
   salinity = T_fabio.PracticalSalinity_psu_;
   density = T_fabio.Density_kg_m3_;
   conduct = T_fabio.Conductivity_S_m_;
   turbid_cell = T_fabio.Turbidity_NTU_;
   cholro = T_fabio.Chlorophyll_ug_l_;
   oxygen = T_fabio.OxygenConcentrationCorrected_ml_l_;

   temp_str = char(temp_cell);
   tempU = temp_str(:,1:5);
   for i=1:length(tempU)
      if ~contains(tempU(i,:),'NA')
         temp(i,:) = str2num(tempU(i,:));
      else
         temp(i,:) = NaN;
      end
   end

   turbid_str = char(turbid_cell);
   turbid_strU = turbid_str(:,1:6);
   for i=1:length(turbid_strU)
      if ~contains(turbid_strU(i,:),'NA')
         turbid(i,:) = str2num(turbid_strU(i,:));
      else
         turbid(i,:) = NaN;
      end
   end
%
%  Define the time
%
   tt = char(timeX);
   tt = tt(:,1:end-1);

   time = datenum(tt,'yyyy-mm-ddTHH:MM:SS');
   time_eq = datenum(2014,09,10,16,36,43);
end % iwant_2014

if iwant_2018_all==1
    T = readtable('fabio_2018.csv');
    loc = 'Fabio: BarkleyCanyon Axis (20160614T090000Z-20191125T185000Z) Avg10minute';

    yy = T.Var1;
    mm = T.Var2;
    dy = T.Var3;
    hr = T.Var4;
    mn = T.Var5;
    sc = T.Var6;

    turbid_orig = T.Var7;

    time_orig = datenum(yy,mm,dy,hr,mn,sc);
   
    %disp('Limiting time to a max of 22-June-2018');
    %ikeep = find(time_orig<datenum(2018,06,22,0,0,0));
    %time = time_orig(ikeep);
    %turbid = turbid_orig(ikeep);

    %disp('Limiting time to a after 22-June-2018');
    %ikeep = find(time_orig>=datenum(2018,06,22,0,0,0));
    %time = time_orig(ikeep);
    %turbid = turbid_orig(ikeep);

    %disp('Limiting time to April - Nov, 2019');
    %ikeepA = find(time_orig<=datenum(2019,11,01,0,0,0));
    %ikeepB = find(time_orig>=datenum(2019,04,15,0,0,0));
    %ikeep = intersect(ikeepA,ikeepB);
    %time = time_orig(ikeep);
    %turbid = turbid_orig(ikeep);

    disp('Keeping all data');
    time = time_orig;
    turbid = turbid_orig;

    tmin = min(time);
    tmax = max(time);
    
end % iwant_2018_all

if iwant_2018_some==1
    T = readtable('fabio_2018.csv');
    loc = 'Fabio: BarkleyCanyon Axis (20160614T090000Z-20191125T185000Z) Avg10minute';

    yy = T.Var1;
    mm = T.Var2;
    dy = T.Var3;
    hr = T.Var4;
    mn = T.Var5;
    sc = T.Var6;

    turbid_orig = T.Var7;

    time_orig = datenum(yy,mm,dy,hr,mn,sc);
   
    %disp('Limiting time to a max of 22-June-2018');
    %ikeep = find(time_orig<datenum(2018,06,22,0,0,0));
    %time = time_orig(ikeep);
    %turbid = turbid_orig(ikeep);

    disp('Limiting time to a after Aug-2018 and before April 2019');
    ikeepA = find(time_orig>=datenum(2018,08,01,0,0,0));
    ikeepB = find(time_orig<=datenum(2019,04,01,0,0,0));
    ikeep = intersect(ikeepA,ikeepB);
    time = time_orig(ikeep);
    turbid = turbid_orig(ikeep);

    %disp('Limiting time to April - Nov, 2019');
    %ikeepA = find(time_orig<=datenum(2019,11,01,0,0,0));
    %ikeepB = find(time_orig>=datenum(2019,04,15,0,0,0));
    %ikeep = intersect(ikeepA,ikeepB);
    %time = time_orig(ikeep);
    %turbid = turbid_orig(ikeep);

    tmin = min(time);
    tmax = max(time);
end
    
%
%
%  Load earthquake data
%
load eqs_over55.dat
eq55_yy =  eqs_over55(:,1); 
eq55_mm =  eqs_over55(:,2);
eq55_dy =  eqs_over55(:,3);
eq55_hr =  eqs_over55(:,4);
eq55_mn =  eqs_over55(:,5);
eq55_sc =  eqs_over55(:,6);
eq55_lat = eqs_over55(:,7);
eq55_lon = eqs_over55(:,8);
eq55_dep = eqs_over55(:,9);
eq55_mag = eqs_over55(:,10);
eq55_time = datenum(eq55_yy,eq55_mm,eq55_dy,eq55_hr,eq55_mn,eq55_sc);
n55 = length(eq55_time);

load eqs_over7.dat
eq7_yy =  eqs_over7(:,1); 
eq7_mm =  eqs_over7(:,2);
eq7_dy =  eqs_over7(:,3);
eq7_hr =  eqs_over7(:,4);
eq7_mn =  eqs_over7(:,5);
eq7_sc =  eqs_over7(:,6);
eq7_lat = eqs_over7(:,7);
eq7_lon = eqs_over7(:,8);
eq7_dep = eqs_over7(:,9);
eq7_mag = eqs_over7(:,10);

eq7_time = datenum(eq7_yy,eq7_mm,eq7_dy,eq7_hr,eq7_mn,eq7_sc);
n7 = length(eq7_time);
%
%  Define favorite data
%
time_fav(1) = datenum(2012, 10, 28,  03,04,08);
time_fav(2) = datenum(2011, 09, 09,  19,41,34);
time_fav(3) = datenum(2018, 10, 22,  06,16,26);
time_fav(4) = datenum(2011, 03, 11,  05,46,24);
time_fav(5) = datenum(2014, 04, 24,  03,10,10);
time_fav(6) = datenum(2013, 01, 05,  08,58,14);
time_fav(7) = datenum(2018, 01, 23,  09,31,40);
time_fav(8) = datenum(2018, 10, 22,  05,39,39);
time_fav(9) = datenum(2005, 06, 15,  02,50,54);
time_fav(10) = datenum(2018, 10, 22,  06,22,48);
time_fav(11) = datenum(2012, 11, 08,  02,01,50);
time_fav(12) = datenum(2008, 01, 05,  11,01,06);
time_fav(13) = datenum(2013, 05, 24,  05,44,48);
time_fav(14) = datenum(2010, 02, 27,  06,34,11);
time_fav(15) = datenum(2017, 09, 08,  04,49,19);
time_fav(16) = datenum(2006, 11, 15,  11,14,13);

nfav = length(time_fav);
disp(['Number of high strain events: ',num2str(nfav)]);
%
%  Loop over all percentages
%
for nnn = 1:length(Pwant);
%
%  Define what percentage to use in the percentile tests
%
   Pat  = Pwant(nnn); 
   disp(['Workingi on ',num2str(Pat),'%']);
%
%  Clear all arrays
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
   
       time_orig = time;
       turbid_orig = turbid;
%
%  Define stuff
%
      N = length(time_orig);

      time_start = datestr(min(time_orig))
      time_end = datestr(max(time_orig))
      disp('Using all data, no need to loop over 3-year time windows ....');
%
% Make sure data is time sorted
%

         [time_sort,isort] = sort(time_orig);
         turbid_sort = turbid_orig(isort);
         N = length(time_sort);
         time = time_sort;
         turid = turbid_sort;
         N = length(time);
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
            while(t2<tend)
                iselect = find(and(time<t2,time>=t1));
                nfound(iat) = sum(turbid(iselect));
                timestamp(iat) = t1;
                %timestamp(iat) = t2;
                %timestamp(iat) = t1+0.5*(t2-t1);
                t1 = t1+tmovit_days;
                t2 = t1+testwin_days;
                iat = iat + 1;
            end
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
timeA = datenum(2018,10,29,18,58,36);
timeB = datenum(2018,11,15,10,32,24);
pp=patch([timeA timeB timeB timeA timeA],[0 0 83.0 83.0 0],'y');
set(pp,'FaceAlpha',0.50,'LineStyle','none');
hold on
               h=plot(timestamp,(nfound),'k.','Color',[0.0 0.6 0.4]);
               hold on
               axis tight
               xx = get(gca,'Xlim');
               plot(timestamp(ioutlie),(nfound(ioutlie)),'ro','Color',[1.0 0.6 0.0],'LineWidth',1)
               %ylabel('Turbidity')
               ylabel('S_T')
               plot(xx,[nthresh nthresh],'--','Color',[1.0 0.6 0.0],'LineWidth',2);
            end
            xlabel('Time')
            yy = get(gca,'Ylim');
            yy = get(gca,'Ylim');
            for i=1:n55
               plot([eq55_time(i) eq55_time(i)],yy,'k:','Color',[0.8 0.8 0.8],'LineWidth',2);
            end
            for i=1:n7
               plot([eq7_time(i) eq7_time(i)],yy,'k--','Color',[0.8 0.8 0.8],'LineWidth',2);
            end
            for i=1:nfav
               plot([time_fav(i) time_fav(i)],yy,'r-','Color',[0.6 0.6 0.4],'LineWidth',3);
            end
 
            %h=text( 737474.7311,      69.53555556, 'Barkley Canyon Axis'); 
            %set(h,'FontName','Helvetica','FontSize',14);
  
            set(gca,'Xlim',[tstart tend]);
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
            title(['Barkley Canyon Axis (',num2str(Pat),'%): \Deltatime=',num2str(testwin_hours),'h, moving shift=',num2str(tmoveit_hours),'h, ',tlab3]);
            tout = [pdir,'/',Pname_base,'_paper.png'];
            tout = [pdir,'/',Pname_base,'_Dt_',num2str(testwin_hours),'_',datestr(min(time),30),'_',datestr(max(time),30),'_Percent_',num2str(Pat),'.png'];
            set(gca,'FontSize',12)
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
