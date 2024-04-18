%
% Process all years of data
%
clear all
close all

P = path;
path(P,'../Matlab_codes_NEHRP');


disp('*** Add in 2019 at some point');
years = [2016, 2017, 2018, 2019, 2020, 2021, 2022];
%years = [2017];
%years = [2016, 2017, 2018, 2019, 2020, 2021, 2022];
%years = [2018, 2019];
%years = [2019];
nyears = length(years)
for i=1:nyears
    yearat = years(i);
    dirat = ['../process_',num2str(yearat),'_new_202401/'];
    if ~exist(dirat,'dir')
        disp('Output directory does not exist, now creating');
        system(['mkdir ',dirat]);
    end
    system(['cd ../',dirat]); % Move to output directory where we want to store everything
    disp(['********** Processing ',num2str(yearat)]);
    sub_process_data_firfilt(yearat,dirat)
    system(['cd ../Matlab_codes_NEHRP']);
end
