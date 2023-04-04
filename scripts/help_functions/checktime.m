function data = checktime(data, config, reader)
% RG 27.3.2023

% 1. check if all timestamps in the correct hour

[yy, mm, dd, hh, ~, ~] = datevec(double(data.time(1))/3600/24 + datenum([2001,1,1,0,0,0])); % assuming that first time is correct
ind = data.time < datetimeconv(yy, mm, dd, hh, 0, 0) | data.time > datetimeconv(yy, mm, dd, hh+1, 0, 0);

test1 = sum(ind) ~= 0;


% 2. check if any dubplicate timestamps  
test2 = any(diff(double(data.time) + double(data.sampleTms).*1e-3) <= 0 );


% 3. if either test is true, continue checking time from lv1 file
if ~(test1 || test2)
    return
end
  

%% reading lv1_file, and comparing time stamp from lv0 and lv1 files


% find lv1 corresponding to lv0 file
lv1_file = config.infile( 1 : strfind(config.infile,config.filetype(2:end))-1); % strip file name (incl. path) without file ending
lv1_file = [lv1_file config.filetype_lv1];
 
% read lv1 file
data_lv1 = reader.lv1(lv1_file);  

if all(data_lv1.time == data.time) % time stamps are identical
    return
end

error('lv1 file timestamps differ from lv0 timestamps.')
