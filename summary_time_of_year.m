T2016 = readtable('2016_plots_filter_1e-06_0.0079167_despike_length_3/Data/Template_list_created_20230126T141427_filtered_2_15_0.18654_1_0_0.2.csv');
T2017 = readtable('2017_plots_filter_1e-06_0.0079167_despike_length_3/Data/Template_list_created_20230126T141440_filtered_2_15_2.0921_1_0_0.2.csv');
T2018 = readtable('2018_plots_filter_1e-06_0.0079167_despike_length_3/Data/Template_list_created_20230126T150249_filtered_2_15_0.20451_1_0_0.2.csv');
T2019 = readtable('2019_plots_filter_1e-06_0.0079167_despike_length_3/Data/Template_list_created_20230126T150636_filtered_2_15_1.136_1_0_0.2.csv');
T2020 = readtable('2020_plots_filter_1e-06_0.0079167_despike_length_3/Data/Template_list_created_20230126T141758_filtered_2_15_0.19284_1_0_0.2.csv');
T2021 = readtable('2021_plots_filter_1e-06_0.0079167_despike_length_3/Data/Template_list_created_20230126T142302_filtered_2_15_1.0099_1_0_0.2.csv');
T2022 = readtable('2022_plots_filter_1e-06_0.0079167_despike_length_3/Data/Template_list_created_20230126T142526_filtered_2_15_1.8382_1_0_0.2.csv');

%time_2016 = T2016.Var1;
time_2017 = T2017.Var1;
time_2018 = T2018.Var1;
time_2019 = T2019.Var1;
time_2020 = T2020.Var1;
time_2021 = T2021.Var1;
time_2022 = T2022.Var1;

figure(1)
clf
%plot(time_2016,16*ones(size(time_2016)),'go');
%hold on
plot(time_2017,17*ones(size(time_2017)),'go');
hold on
plot(time_2018,18*ones(size(time_2018)),'go');
plot(time_2019,19*ones(size(time_2019)),'go');
plot(time_2020,20*ones(size(time_2020)),'go');
plot(time_2021,21*ones(size(time_2021)),'go');
plot(time_2022,22*ones(size(time_2022)),'go');
