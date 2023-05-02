function data = checktime(data, config, reader)
% RG 27.3.2023

% 1. check if all timestamps in the correct hour

[yy, mm, dd, hh, ~, ~] = datevec(double(data.time(1))/3600/24 + datenum([2001,1,1,0,0,0])); % assuming that first time is correct
ind = data.time < datetimeconv(yy, mm, dd, hh, 0, 0) | data.time > datetimeconv(yy, mm, dd, hh+1, 0, 0);

test1 = sum(ind) ~= 0;


% 2. check if any dubplicate timestamps 
ind2 = diff(double(data.time) + double(data.sampleTms).*1e-3) <= 0 ;
test2 = any(ind2);


% 3. if either test is true, continue checking time from lv1 file
if ~(test1 || test2)
    return
end
  

%% reading lv1_file, and comparing time stamp from lv0 and lv1 files


% find lv1 corresponding to lv0 file
lv1_file = config.infile( 1 : strfind(config.infile,config.filetype(2:end))-1); % strip file name (incl. path) without file ending
lv1_file = [lv1_file config.filetype_lv1];

% check if lv1 file reader available
if ~isfield(reader, 'lv1')
    return
end

% check if lv1 file exists
if ~exist(lv1_file, 'file')
    return
end

% read lv1 file
data_lv1 = reader.lv1(lv1_file);  

if isempty(data_lv1)
    return
end

% check if time arrays the same size, if not, nothing to be done
if length(data_lv1.time) ~= length(data.time)
    return
end

if all(data_lv1.time == data.time) % time stamps are identical
    return
end

% checked, lv1 file stamps does not look better than lv0
%if strcmp( lv1_file, '/data/obs/site/nya/joyrad94/l0/2022/04/26/joyrad94_20220426040003_P01_ZEN.lv1')
%    return
%end

% check if time in lv1 file has problems
ind_test1 = data_lv1.time < datetimeconv(yy, mm, dd, hh, 0, 0) | data_lv1.time > datetimeconv(yy, mm, dd, hh+1, 0, 0);
ind_test2 = diff(double(data_lv1.time) + double(data_lv1.sampleTms).*1e-3) <= 0 ;

% if no problems in lv 1 file, then it is definitely better
if sum(ind_test1) == 0 && sum(ind_test2) == 0
    data.time = data_lv1.time;
    data.sampleTms = data_lv1.sampleTms;
    return

% check if less problematic time stamps in lv1 file than lv0 file
elseif ( sum(ind_test1) + sum(ind_test2) ) < ( sum(ind) + sum(ind2) )
    data.time = data_lv1.time;
    data.sampleTms = data_lv1.sampleTms;
    return
    
else 
    disp('checked time stamp from l1 file, not better than in the lv0 file -> not considering them further.')
    
end

