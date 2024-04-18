%
% Goal: Take a look at the base data.
% Motivation: Filtering produced questionable results, need to figure out why
%
clear all
close all

P = path;
path(P,'../Matlab_codes_NEHRP');


years = [2016, 2017, 2018, 2019, 2020, 2021, 2022];

nyears = length(years)
disp(['Looping over N=',num2str(nyears),' years of data.']);
for i=1:nyears
    yearat = years(i);
    dirat = ['../data_check_202401/'];
    if ~exist(dirat,'dir')
        disp('Output directory does not exist, now creating');
        system(['mkdir ',dirat]);
    end
    system(['cd ../',dirat]); % Move to output directory where we want to store everything
    disp(['********** Processing ',num2str(yearat)]);
    sub_see_data_firfilt(yearat,dirat)
    system(['cd ../Matlab_codes_NEHRP']);
end
