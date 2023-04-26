function data = checkdata1(data)
% sometimes data logging problems - remove these from the data
% using time stamp for check, if time is outside the hour for hour of the
% data, remove these instances from the all time series
% RG 14.8.2019        


[yy, mm, dd, hh, ~, ~] = datevec(double(data.time(1))/3600/24 + datenum([2001,1,1,0,0,0])); % assuming that first time is correct

ind = data.time < datetimeconv(yy, mm, dd, hh, 0, 0) | data.time > datetimeconv(yy, mm, dd, hh+1, 0, 0);

iout = find(ind);

if isempty(iout)
    return
end


dim_t = length(data.time);

% loop for all variables
fields = fieldnames(data);

for ff = 1:length(fields)

    vardim = size(data.(fields{ff}));

    if ~any(vardim == dim_t)
        continue
    end


    if strcmp(fields{ff}, 'QualFlag')
        data.(fields{ff})(iout,:,:) = [];

    elseif prod(vardim) == dim_t % 1 dim variables
        data.(fields{ff})(iout) = [];

    elseif length(vardim) == 2 % 2 dim variables
        data.(fields{ff})(iout,:) = [];

    elseif length(vardim) == 3 % 3 dim variables
        data.(fields{ff})(iout,:,:) = [];

    end
end


data.totsamp = length(data.time); 
data.totsampchangelabel = 1; % note for later that issues with time found
