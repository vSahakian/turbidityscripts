function sub_process(yearat,dirat)
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
Mname_filt = [dirat,num2str(yearat),'_metadata_filtered.csv']
fidM_filt = fopen(Mname_filt,'w+');

Mname_nofilt = [dirat,num2str(yearat),'_metadata_unfiltered.csv']
fidM_unfilt = fopen(Mname_nofilt,'w+');
%
%  Define subdirectories for plots and data
%
Oname_data = [dirat,'data_',num2str(yearat),'/'];
Oname_plots = [dirat,'plots_',num2str(yearat),'/'];
system(['mkdir ',Oname_plots]);
system(['mkdir ',Oname_data]);
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
disp(['Now processing ',unique(y),' data.']);
%
%  Remove spikes using medfilt1
%
despike_len = 3;
data_despike = medfilt1(data,despike_len);
%
%  Labels for data we will toss
%
toss_lab = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'];
%
%  Plot oringal data and despiked data
%
figure(1)
clf
subplot(3,1,1)
plot(data)
axis tight
title('Original Data');
ylabel('Turbidity');
subplot(3,1,2)
plot(data_despike)
axis tight
ylabel('Turbidity');
title(['Data despiked, len=',num2str(despike_len)]);
subplot(3,1,3)
plot(data-data_despike);
axis tight
xlabel('Index');
ylabel('Residule');
%
%  Filter data
%
disp('Filtering!');
filt_low = 1e-5;
filt_high = 0.95*1/(2*60); % set to nyquest
%OLD DO NOT USE data_filt = eqfiltfilt(data_despike,filt_low,filt_high,delta_sec,4);
data_filt = eqbandpassfir(data_despike,filt_low,filt_high,delta_sec,4);
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
%
%  Compute moving mean in order to establish the noise floor
%
MMweek = movmean(data_filt,tstep_week); % these are arrays
SSweek = movstd(data_filt,tstep_week);
TThresh_week = MMweek + 2*SSweek;
%
% Plot raw data
%
figure(2)
clf
plot(time_days,data,'k-','Color',[0.8 0.8 0.8]);
hold on
plot(time_days,data_filt,'g-','Color',[0.2 0.8 0.4]);
hold on
axis tight
title(['Data ',num2str(yearat),': tstep=',num2str(round(tstep_week_days)),'days; filtbp: ',num2str(filt_low),'-',num2str(filt_high)])
grid on
ylabel('Turbidity (units)')
xlabel('Time')
set(gca,'FontSize',12)
yy = get(gca,'Ylim');
plot(time_days,MMweek,'k-','Color',[0.8 0.6 0.4],'LineWidth',2);
datetick('x');
legend('Original Data','Filtered Data',['Running mean of consecutive points (N=',num2str(round(tstep_week_days)),' days)'],'Location','southoutside');
set(gca,'FontSize',11);
pname = [dirat,Oname_plots,'plot_data_',num2str(yearat),'.png'];
print(gcf,pname,'-dpng');
%
%  Define STA / LTA parameters
% 1. Short term length
% 2. Long term length
% 3. Threshold to trigger
% 4. Trigger off when ratio drops below this value
%
%  Process each tstep duration of data separately using time in minutes
%
i1 = 1;
istop = 0;
tdetect_min = [];
idetect1_days = [];
idetect2_days = [];
while istop==0
    i2 = i1+tstep_week;
    if i2>=length(data_filt)
        i2 = length(data_filt);
        istop = 1;
    end
    iindex = [i1:i2];
    %edp = [1,3, TThresh_week(i1),1.0,0,0.5];
    %edp = [2,10, TThresh_week(i1),1.0,0,0.5];
    %edp = [2,7, TThresh_week(i1),1.0,0,0.5];
    %edp = [2,15, TThresh_week(i1),1.0,0,0.5];
    %edp = [2,15, TThresh_week(i1),1.0,0,0.25];
    edp = [2,15, TThresh_week(i1),1.0,0,0.20];
    edp_lab = [num2str(edp(1)),'_',num2str(edp(2)),'_',num2str(edp(3)),'_',num2str(edp(4)),'_',num2str(edp(5)),'_',num2str(edp(6))];
    %
    % define time snippet
    %
    time_rel_min_snip = time_rel_min(i1:i2);
    data_snip = data_filt(i1:i2);
    %
    % Is there a signal in this time snippit?
    %
    tdetect_temp = sta_lta_2(time_rel_min_snip,data_snip,edp,'sst'); % sst requests the start and stop time in mins
    idetect_temp = sta_lta_2(time_rel_min_snip,data_snip,edp,'ssd'); % ssd returns the start and stop index
    %
    % Save time and index information if a signal was identified
    %
    if length(tdetect_temp)>0
        tdetect_min = [tdetect_min;tdetect_temp(:,1)]; % actual time in min
    end
    if length(idetect_temp)>0
        idetect1_days = [idetect1_days,iindex(idetect_temp(:,1))];
        idetect2_days = [idetect2_days,iindex(idetect_temp(:,2))];
    end
    i1 = i1+tstep_week;
end
disp(['Number of detected signals=',num2str(length(idetect1_days))]);

figure(3)
clf
plot(time_days,data_filt,'g-','Color',[0.8 0.6 1.0]) % Replot so above trigger marks
hold on
axis tight
yy = get(gca,'Ylim');
for i=1:length(idetect1_days)
    plot([time_days(idetect1_days(i)) time_days(idetect1_days(i))],yy,'k-','Color',[0.4 0.6 0.8]);
end
grid on
edpout = strrep(num2str(edp),'      ',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
title({['Number of Triggers = ',num2str(length(idetect1_days)),' (parmameters=',edpout,')'],t2});
datetick('x');
xlabel('Time')
ylabel('Turbidity (units)');
set(gca,'FontSize',12)
pname = [dirat,Oname_plots,'plot_triggers',num2str(yearat),'.png'];
print(gcf,pname,'-dpng');
%
%  Loop over all detections to see if they have the proper fall off
%
igood = 0;
ibad = 0;
figure(4)
clf
icount = 0;
ntossA = 0;
indexA = [];
disp('Starting loop');
iflag_save = NaN;
iadd = 100; % Number of data points to investigate
numdup = 0;
itemplate = 0;
for i=1:length(idetect1_days);
    flag(i) = 0; % set to zero to start ....
    pdup = '';
    i1=find(time_rel_min>=tdetect_min(i,1)); % List of index at or above dection
    i2 = find(time_rel_min<tdetect_min(i,1)+iadd); % list of times below tdetect + iadd data points
    iwant = intersect(i1,i2);

    time_snip_days_orig = time_days(iwant);
    time_snip_orig = time_rel_min(iwant);
    data_snip_orig = data_filt(iwant);
    %
    % Add iadd points prior to the detection
    %
    ipre = max(min(iwant)-iadd,1);
    time_snip_longer = time_rel_min(ipre:max(iwant));
    time_snip_longer_days = time_days(ipre:max(iwant));
    data_snip_longer = data_filt(ipre:max(iwant));
    data_snip_longer_nofilt = data(ipre:max(iwant));
    %
    %  Determime the time of the maximum signal
    %
    [ii,jj] = max(data_snip_orig);
    tsave(i) = time_snip_orig(jj);
    tsave_days(i) = time_snip_days_orig(jj);

    time_snip = time_snip_orig-time_snip_orig(jj); % shift to have maximum at time = 0;
    data_snip = data_snip_orig;

    time_snip_longer = time_snip_longer-time_snip_orig(jj); % shift longer waveform to 0 too

    [iii,jjj] = max(data_snip);
    if jjj>=40
        % disp('Max can not be the end');
        %[iii,jjj] = max(data_snip(1:25));
        % The above forces the length of the data snip < 40 so will be flaged below
        ibigend(i) = 1;
    else
        ibigend(i) = 0;
    end

    x = time_snip(jjj:end);
    y = data_snip(jjj:end);
    %
    %  Require a certain number of points in the signal to be over 2 standard deviations from the mean
    %
    thresh = MMweek(jjj) + 2*SSweek(jjj);
    itemp = find(y>thresh);
    iover = length(itemp);
    if iover>=25
        iambig(i) = 1;
    else
        iambig(i) = 0;
    end
    %
    %   Check if this signal was already detected
    %
    if isnan(iflag_save)
        inotdup = 1; % First detection is imposible to be a dup
    elseif ((min(iwant)+jjj)-iflag_save) < 50 % duplicate
        %xx = get(gca,'Xlim');
        %yy = get(gca,'Ylim');
        %text(xx(1) + 0.2*(xx(2)-xx(1)),yy(2)-0.4*(abs(yy(1)-yy(2))),'Duplicate?')
        pdup = 'duplicate_';
        inotdup = 0;
        iflag_save = min(iwant)+jjj; % Reset iflag_save anyway .... ???? right ???
        numdup = numdup+1;
    else % not duplicate
        inotdup = 1;
    end
    %
    %  Check if data good enought to fit
    %
    if and(inotdup,and(ibigend(i)==0,iambig(i)==1))
        %disp('DEBUG: data good enough to fit');
        f0 = fit(x,y,'exp1','StartPoint',[max(y),-.1]); % exp1 ==> Y = a*exp(b*x)
        c95 = confint(f0);
        if and(-0.2<=f0.b,f0.b<=-0.02);
            icount = icount + 1;
            if ifitfig == 1
                figure(4)
                subplot(2,1,1)
                plot([icount icount],[(c95(1,1)) (c95(2,1))],'r-','LineWidth',1);
                hold on
                plot(icount,f0.a,'rx','MarkerSize',15,'LineWidth',1,'Color',[0.7 0.4 0.2]);
                subplot(2,1,2)
                plot([icount icount],[c95(1,2) c95(2,2)],'b-','LineWidth',1);
                hold on
                plot(icount,f0.b,'bo','MarkerSize',15,'LineWidth',1);
            end % check ifull
            %
            %  Check that the signal identified larger
            %
            ipre = intersect(find(time_snip_longer>-30),find(time_snip_longer<0));
            ipost = find(time_snip_longer>0);
            mpre = median(data_snip_longer(ipre));
            mpost = median(data_snip_longer(ipost));
            if mpre>mpost
                %disp('Nope!')
                itoss(i) = 1;
            else
                itoss(i) = 0;
            end
            if itoss(i)==0
                %disp('We have a contender, lets check it out!');
                figure(5)
                clf
                l1=plot(time_snip_longer,data_snip_longer,'ko--');
                hold on
                l1X = plot(time_snip_longer,data_snip_longer_nofilt,'ko-','Color',[0.8 0.8 0.8]);
                ytopA = max(data_snip_longer) + 0.1*max(data_snip_longer);
                ytopB = max(data_snip_longer_nofilt) + 0.1*max(data_snip_longer_nofilt);
                ytop = max(ytopA,ytopB);
                l2=plot(x,y,'bo-','Color',[0.0 0.8 0.4]);
                hold on
                xx = time_snip;
                l3=plot(xx,f0(xx),'r-');
                labdelta = ['\DeltaB',num2str(abs(c95(1,2)-c95(2,2)),2),' \DeltaA=',num2str(abs(c95(1,1)-c95(2,1)),2)];
                title({['FLAG #',num2str(icount),' ',datestr(tsave_days(i),30),' b=',num2str(f0.b,'%3.2f'),'(',num2str(c95(1,2),2),'-',num2str(c95(2,2),2),') ',labdelta,' Thresh ',num2str(thresh,2),' iover=',num2str(iover),' iambig=',num2str(iambig(i))],t2});
                amp_sig(icount) = max(y);
                time_sig(icount) = tsave_days(i);
                thresh_sig(icount) = thresh;
                flag(i) = 1;
                grid on
                axis tight
                xx = get(gca,'Xlim');
                l4=plot(xx,[thresh thresh],'k--');
                legend([l1,l1X,l2,l3,l4],'Filtered data','Data',['Data, N=',num2str(length(x))],'Fit Curve','mean+2*std');
                set(gca,'Ylim',[0 ytop]);
                set(gca,'FontSize',10);
                if abs(c95(1,1)-c95(2,1))>8
                    ntossA = ntossA + 1;
                    indexA(ntossA) = i;
                    ismember(i,indexA)
                    xx = get(gca,'Xlim');
                    yy = get(gca,'Ylim');
                    text(xx(1) + 0.2*(xx(2)-xx(1)),yy(2)-0.2*(abs(yy(1)-yy(2))),'Toss this because \DeltaA>8?')
                end
                if itoss(i)==1
                    xx = get(gca,'Xlim');
                    yy = get(gca,'Ylim');
                    text(xx(1) + 0.2*(xx(2)-xx(1)),yy(2)-0.4*(abs(yy(1)-yy(2))),'Toss because in the shadow of prior signal?')
                end
                %
                %     Save index to check for dups on next loop
                %
                iflag_save = min(iwant) + jjj; % save this index so we can avoid duplicates
                iflag(i) = 1; % flag to indicate this is a good detection
                %
                %  Ask to save
                %
                %ans = input('Keep [Y/N]?','s');
                disp('Hard coded!');
                ans = 'y';
                if or(ans=='Y',ans=='y');
                    %disp('Bingo -- keeping');
                    itemplate = itemplate+1;
                    figure(5)
                    xx = get(gca,'Xlim');
                    yy = get(gca,'Ylim');
                    text(xx(1) + 0.2*(xx(2)-xx(1)),yy(2)-0.4*(abs(yy(1)-yy(2))),['Save as template #',num2str(itemplate)])
                    time_template(itemplate) = time_sig(icount); % save favorite times
                    %
                    % Create csv files for filtered data, including metadata
                    %
                    %fname = [Oname_data,'template_filtered_',num2str(filt_low),'_',num2str(filt_high),'_',edp_lab,'_',datestr(time_sig(icount),30),'.csv'];
                    %fname = [dirat,Oname_data,'template_',num2str(icount,'%03.f'),'_filtered_',datestr(time_sig(icount),30),'.csv'];
                    fname = [dirat,Oname_data,'template_',num2str(itemplate,'%03.f'),'_filtered_',datestr(time_sig(icount),30),'.csv'];
                    fid = fopen(fname,'w+');
                    for iii=1:length(time_snip_longer_days)
                        tout = [datestr(time_snip_longer_days(iii),30),', ',num2str(data_snip_longer(iii))];
                        fprintf(fid,'%s\n',tout);
                    end
                    fclose(fid);
                    %
                    if itemplate==1
                        fprintf(fidM_filt,'%s\n',['#   year: ',num2str(yearat)]);
                        fprintf(fidM_filt,'%s\n',['#   dt: ',num2str(dt)]);
                        fprintf(fidM_filt,'%s\n',['#   despike_len: ',num2str(despike_len)]);
                        fprintf(fidM_filt,'%s\n',['#   number of days used to establish noise floors: ',num2str(tstep_week_days)]);
                        fprintf(fidM_filt,'%s\n',['#   lta_sta_params: ',edp_lab]);
                        fprintf(fidM_filt,'%s\n',['#   filter band: ',num2str(filt_low),'-',num2str(filt_high)]);
                        fprintf(fidM_filt,'%s\n',['template#,time_min,time_max,time_event,amp_event,fname']);
                    end
                    mout = ['template_',num2str(itemplate,'%03.f'),',',datestr(time_snip_longer_days(1),30),',',datestr(time_snip_longer_days(end),30),',',datestr(time_sig(icount),30),',',num2str(amp_sig(icount)),',',fname];
                    fprintf(fidM_filt,'%s\n',mout);
                    %
                    % Create csv files for non-filtered data, including metadata
                    %
                    fname = [dirat,Oname_data,'template_',num2str(itemplate,'%03.f'),'_notfiltered_',datestr(time_sig(icount),30),'.csv'];
                    fid = fopen(fname,'w+');
                    for iii=1:length(time_snip_longer_days)
                        tout = [datestr(time_snip_longer_days(iii),30),', ',num2str(data_snip_longer_nofilt(iii))];
                        fprintf(fid,'%s\n',tout);
                    end
                    fclose(fid);
                    %
                    % Now add metadata file
                    %    name,StartTime,EventStartTime,EndTime,path,filt_low,filt_high,despikeLen,edp_lab,RunDate,
                    %time_days(min(iwant))
                    %time_days(max(iwant));
                    %
                    if itemplate==1
                        fprintf(fidM_unfilt,'%s\n',['#   year: ',num2str(yearat)]);
                        fprintf(fidM_unfilt,'%s\n',['#   dt: ',num2str(despike_len)]);
                        fprintf(fidM_unfilt,'%s\n',['#   despike_len: ',num2str(despike_len)]);
                        fprintf(fidM_unfilt,'%s\n',['#   number of days used to establish noise floors: ',num2str(tstep_week_days)]);
                        fprintf(fidM_unfilt,'%s\n',['#   lta_sta_params: ',edp_lab]);
                        fprintf(fidM_unfilt,'%s\n',['#   {Not filtered}']);
                        fprintf(fidM_unfilt,'%s\n',['template#,time_min,time_max,time_event,amp_event,fname']);
                    end
                    mout = ['template_',num2str(itemplate,'%03.f'),',',datestr(time_snip_longer_days(1),30),',',datestr(time_snip_longer_days(end),30),',',datestr(time_sig(icount),30),',',num2str(amp_sig(icount)),',',fname];
                    fprintf(fidM_unfilt,'%s\n',mout);
                    ireview(icount) = 1;
                else
                    itoss(i) = 2; % tossing after manual review
                    ireview(icount) = 0;
                    disp('Tossing');
                end
            else
                iplotbad = 0;
                if iplotbad
                    figure(5)
                    clf
                    l1=plot(time_snip_longer,data_snip_longer,'ko-','Color',[0.8 0.8 0.8]);
                    hold on
                    l2=plot(x,y,'ko-');
                    hold on
                    l3=plot(xx,f0(xx),'r-');
                    title(['Detection #',num2str(i),'@time=',num2str(tsave(i)),' b=',num2str(f0.b,'%3.2f'),' not flagged; iambig=',num2str(iambig(i))']);
                    ibad = ibad + 1;
                    grid on
                    xx = get(gca,'Xlim');
                    plot(xx,[thresh thresh],'k--');
                    set(gca,'FontSize',12);
                end
                flag(i) = 0;
            end
            %legend([l1,l2,l3],'Extended data',['Data, N=',num2str(length(xx))],'Fit Curve');
            figure(5)
            grid on
            xlabel('Time (min)');
            ylabel('Turbidity (UNITS)');
            if flag(i) == 1
                if find(ismember(i,indexA))
                    pname = [dirat,Oname_plots,'plot_',num2str(yearat),'_template_',toss_lab(ntossA),'_tossA.png'];
                else
                    pname = [dirat,Oname_plots,'plot_',num2str(yearat),'_template_',pdup,num2str(itemplate,'%03.f'),'.png'];
                end
            end
            print(gcf,pname,'-dpng');
        end
    end % end loop check on iambig and ibigend
end % end loop over i
disp('End loop over i');
%
%  Close all open files
%
fclose(fidM_filt);
fclose(fidM_unfilt);
%
%  Define the winners
%
disp(['************ #good auto found= ',num2str(icount),'**** igood for sure=',num2str(itemplate),' ***** #bad = ',num2str(ibad),' **** #dup=',num2str(numdup)]);
ipick = find(flag==1);
%
%  Now add where approved flags are ...
%
if ifitfig==1
    figure(4)
    subplot(2,1,1)
    grid on
    title('Red=a & 95%');
    xlabel('Flag Index');
    ylabel('Fit parameter with error bars');
    set(gca,'FontSize',12)
    subplot(2,1,2)
    grid on
    title('Blue=b & 95%');
    xlabel('Flag index');
    ylabel('Fit parameter with error bars');
    set(gca,'FontSize',12)
    pname = [dirat,Oname_plots,'plot_fit_bars_',num2str(yearat),'.png'];
    print(gcf,pname,'-dpng')
end
%
%  Idetify time periods with the noise floor was too high to detect anythning
%
delta_trig = diff(tdetect_min);
ihighnoise=find(delta_trig>1e4);
istart_blindspot = ihighnoise;
iend_blindspot = ihighnoise+1;

%
% Now replot with just the requested triggers
%
figure(7)
clf
plot(time_days,data_filt,'g-','Color',[0.8 0.6 1.0]) % Replot so above trigger marks
hold on
axis tight
yy = get(gca,'Ylim');
for i=1:length(ipick)
    plot([tsave_days(ipick(i)) tsave_days(ipick(i))],yy,'r-','LineWidth',2,'Color',[1.0 0.8 0.6]);
end
grid on
edpout = strrep(num2str(edp),'      ',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
edpout = strrep(edpout,',,',',');
title({['Number of Triggers approved = ',num2str(itemplate),' (parmameters=',edpout,')'],t2});
datetick('x');
xlabel('Time')
ylabel('Turbidity (units)');
set(gca,'FontSize',12)
yy=get(gca,'Ylim');
pname = [dirat,Oname_plots,'plot_trigggers_approved_',num2str(yearat),'.png'];
for i=1:length(istart_blindspot);
    tt1 = time_days(istart_blindspot(i));
    tt2 = time_days(iend_blindspot(i));
    h=patch([tt1 tt2 tt2 tt1 tt1],[yy(1) yy(1) yy(2) yy(2) yy(1)],'y');
    set(h,'FaceAlpha',0.25);
end
print(gcf,pname,'-dpng');
%
% Now create some examples of B
%   ans(a,b,x) = a*exp(b*x)
%
figure(8)
clf
b = -1*[0.001 0.002 0.003 0.004 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1.0 2  10.0 ];
x = [0:1:100];
a = 100;
cc = linspecer(length(b));
for i=1:length(b)
    f = a*exp(b(i).*x);
    if or(i<=6,i>10)
        plot(x,f,'k--','Color',cc(i,:),'LineWidth',2)
    else
        plot(x,f,'k-','Color',cc(i,:),'LineWidth',4)
    end
    lab(i,:) = ['b=',num2str(b(i),'%07.3f')];
    hold on
end
xlabel('Time')
ylabel('Turbidity')
title('We assume the fall off needs to be similar to the solid traces')
legend(lab)
pname = [dirat,Oname_plots,'plot_bark_fall_off.png'];
print(gcf,pname,'-dpng');

figure(9)
clf
disp('**** assumes 1 second data steps');
ntimes = idetect2_days - idetect1_days; % These are indexes
histogram(ntimes,[min(ntimes):1:100])
mm = mean(ntimes);
ss = std(ntimes);
hold on
axis tight
grid on
yy = get(gca,'Ylim');
l0=plot([mm mm],yy,'r-','LineWidth',2);
l1=plot([(mm+ss) (mm+ss)],yy,'r:','LineWIdth',2);
l2=plot([(mm+2*ss) (mm+2*ss)],yy,'r--','LineWIdth',2);
temp = ['mean (',num2str(mm,2),' min)'];
legend([l0,l1,l2],temp,'mean+1std','mean+2std','Location','NorthEast');
xlabel('Duration of sta/lta detections (min)');
ylabel('Number');
title([num2str(yearat),': Values will depend on sta/lta parameter values'])
set(gca,'FontSize',12);
pname = [dirat,Oname_plots,'detected_signal_durations_',num2str(yearat),'.png'];
print(gcf,pname,'-dpng');
%
%  Save results for later use
%
%save savefor_AGU.mat

%
%  Plot all templates
%
day2min = 1440;
min2day = 1/day2min;
tpre = 30*min2day; % 20 mins prior to peak
tpost = 40*min2day;
figure(10)
clf
pnum = 1;
iat = 1;
for i=1:itemplate
    if iat>5*4
        pname = [dirat,Oname_plots,'templates_',num2str(yearat),'_page',num2str(pnum),'.png'];
        print(gcf,pname,'-dpng');
        figure
        clf
        pnum = pnum + 1;
        iat =  1;
    end
    subplot(5,4,iat)
    iat = iat + 1;
    i1 = find(time_days>=time_template(i)-tpre); % time in days
    i2 = find(time_days<=time_template(i)+tpost);
    iwant = intersect(i1,i2);
    plot(time_days(iwant),data(iwant),'k-','Color',[0.8 0.8 0.8]);
    hold on
    plot(time_days(iwant),data_filt(iwant),'g-');
    grid on
    title(['#',num2str(i),': ',datestr(time_template(i))]);
    axis tight
    set(gca,'FontSize',8);
    xx = get(gca,'Xlim');
    set(gca,'xtick',xx);
    datetick('x','KeepLimits');
end
pname = [dirat,Oname_plots,'templates_',num2str(yearat),'_page',num2str(pnum),'.png'];
print(gcf,pname,'-dpng');
%
%  Now plot the sepctra
%
disp(['DEBUG:',Oname_plots])
sub_spectra(time_days,data,data_filt,Oname_plots);
