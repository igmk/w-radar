function data = checkdata2(data)
% Sometimes issues with time stamp due to GPS failures. Checking if time
% goes "backwards" and correcting these time stamps by linearly
% interpolating between neighbouring times.
% RG & LP 16.8.2019        

if ~any(diff(data.time) <= 0)
    return
end

data.time = double(data.time); % in case time is still an int -> the rest of this function would not work

ind = find(diff(data.time) <= 0);
data.time(ind+1) = NaN;

indnan = isnan(data.time);

xarray = 1:length(data.time);

data.time(indnan) = round(interp1(xarray(~indnan), double(data.time(~indnan)), xarray(indnan) ));

if isnan(data.time(end)) % if last time stamp is same as previous, it does not get interpolated by the function above - replacing none with value before so that next block deals with this
    data.time(end) = data.time(end-1);
end
    
% data.time should be integers, but it is possible that interpolation
% creates real values that are rounded to integers that already exist ->
% non-unique time stamps
Ncnt = 0;

while length(unique( data.time)) ~=  length(data.time)
    Ncnt = Ncnt + 1;

    ind = find(diff(data.time) <= 0);

    lentime = length(data.time);

    % if duplicate or decreasing time is found, removing the data points that occur first in the file
    % added by RG 19.11.2019
    if ~isempty(ind)

        fn = fieldnames(data);
        for k=1:numel(fn)

            if any(size(data.(fn{k})) == lentime)

                timdim = find(size(data.(fn{k})) == lentime);

                switch length(size(data.(fn{k})))

                    case 1
                        data.(fn{k})(ind) = [];

                    case 2
                         switch timdim
                             case 1
                                 data.(fn{k})(ind,:) = [];

                             case 2
                                 data.(fn{k})(:,ind) = [];
                         end

                    case 3

                        switch timdim
                            case 1
                                data.(fn{k})(ind,:,:) = [];
                            case 2
                                data.(fn{k})(:,ind,:) = [];
                            case 3
                                data.(fn{k})(:,:,ind) = [];
                        end
                end
           

            end

        end
    end 


    if Ncnt > 10
       disp('exciting while loop after 10 interations in function checkdata2 - check why issue not resolved')
       break
    end
end



if length(data.time) ~= data.totsamp % change in time dimension 
    data.totsamp = length(data.time);
    data.totsampchangelabel = 1;

end
