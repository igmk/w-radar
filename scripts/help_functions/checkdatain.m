function data = checkdatain(data, config, reader)


% testing if lv1 file can help with the time stamp issue
data = checktime(data, config, reader);

% sometimes data logging problems - remove these from the data
% using time stamp for check, if time is outside the hour for hour of the
% data, remove these instances from the all time series
% RG 14.8.2019        
data = checkdata1(data);


% Sometimes issues with time stamp due to GPS failures. Checking if time
% goes "backwards" or duplicate time stamps presents
% RG & LP 16.8.2019        
data = checkdata2_v2(data);


% Sometimes when measurements are interrupted, the last profile of the 
% file contains nonphysical values - remove all profiles where 
% Ze > 50 dB is found
% RG 7.9.2021        
data = zesanitycheck(data);


% found some cases where status got unreasonable values, these happened in 
% files where started to use lv1 file timestamp, so problaly bad data
% logged and not cathed before
data = statuscheck(data);

% Sometimes in the temperature variables found 0 values, replacing these
% with NaNs since a temperature of 0 K is obviously wrong
% RG 20.10.2022
data = temperaturecheck(data);
