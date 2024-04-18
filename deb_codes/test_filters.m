%
% Goal: Explore different filtering options for the 2023 data
% base_look_all_years.m
%   In particular examine how the filters work for time snippit
%     '02-Jan-2022 14:37:32'
%    '08-Jan-2022 18:02:15'
%
clear all
close all

P = path;
path(P,'../Matlab_codes_NEHRP');


%years = [2016, 2017, 2018, 2019, 2020, 2021, 2022];
years = [2022];

nyears = length(years)
disp(['Looping over N=',num2str(nyears),' years of data.']);
dirat = ['../test_filters/'];
if ~exist(dirat,'dir')
   disp('Output directory does not exist, now creating');
   system(['mkdir ',dirat]);
end
for i=1:nyears
    yearat = years(i);
    system(['cd ../',dirat]); % Move to output directory where we want to store everything
    disp(['********** Processing ',num2str(yearat)]);
    %sub_see_data(yearat,dirat)
    sub_test_filters(yearat,dirat)
    system(['cd ../Matlab_codes_NEHRP']);
end
