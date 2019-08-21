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

data.time(indnan) = round(interp1(xarray(~indnan), data.time(~indnan), xarray(indnan) ));

    
% data.time should be integers, but it is possible that interpolation
% creates real values that are rounded to integers that already exist ->
% non-unique time stamps

while length(unique( data.time)) ~=  length(data.time)
    ind = find(diff(data.time) == 0);
    data.time(ind) = data.time(ind) -1;
end
