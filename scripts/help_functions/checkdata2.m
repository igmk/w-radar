function data = checkdata2(data)
% Sometimes issues with time stamp due to GPS failures. Checking if time
% goes "backwards" and correcting these time stamps by linearly
% interpolating between neighbouring times.
% RG & LP 16.8.2019        

if ~any(diff(data.time) <= 0)
    return
end

ind = find(diff(data.time) <= 0);
data.time(ind+1) = NaN;

indnan = isnan(data.time);

xarray = 1:length(data.time);

data.time(indnan) = interp1(xarray(~indnan), data.time(~indnan), xarray(indnan) );

    