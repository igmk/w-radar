function data = checkdata2_v2(data)
% Sometimes issues with time stamp due to GPS failures. Checking if time
% goes "backwards" or duplicate time stamps presents
% RG 21.10.2022     

if ~any(diff(double(data.time) + double(data.sampleTms).*1e-3) <= 0 ) 
    return
end

timestamp = double(data.time) + double(data.sampleTms).*1e-3;

%% 1st check: is there an isolated forward jumping time stamp?
% remove this
ind = find(diff(timestamp) < 0 );

for tt = fliplr(ind)
    
    % check if next time stamp smaller AND if next time stamp larger than previous one 
    if timestamp(tt+1) < timestamp(tt) && timestamp(tt+1) > timestamp(tt-1)
        % remove all data related to this time stamp
        timestamp(tt) = [];
        data = remove_data(data, tt);
    end
            
end % tt

if length(data.time) ~= data.totsamp % change in time dimension 
    data.totsamp = length(data.time);
    data.totsampchangelabel = 1;

end

% if problem solved, return
if ~any(diff(double(data.time) + double(data.sampleTms).*1e-3) <= 0 ) 
    return
end


%% 2nd check: remove duplicate time stamps

ind = find(diff(timestamp) == 0 );

% remove second instance of the duplicate time stamp
ind = ind + 1 ;
data = remove_data(data, ind);


if length(data.time) ~= data.totsamp % change in time dimension 
    data.totsamp = length(data.time);
    data.totsampchangelabel = 1;

end

% if problem solved, return
if ~any(diff(double(data.time) + double(data.sampleTms).*1e-3) <= 0 ) 
    return
end


%% 3rd check: time jumps backwards
disp('WARNING ! This not programmed yet')


end % function



function data = remove_data(data, ix)

    dim_t = length(data.time);
    
    % loop for all variables
    fields = fieldnames(data);

    for ff = 1:length(fields)

        vardim = size(data.(fields{ff}));

        if ~any(vardim == dim_t)
            continue
        end


        if prod(vardim) == dim_t % 1 dim variables
            data.(fields{ff})(ix) = [];

        elseif length(vardim) == 2 % 2 dim variables
            data.(fields{ff})(ix,:) = [];

        elseif length(vardim) == 3 % 3 dim variables
            data.(fields{ff})(ix,:,:) = [];

        end
    end
    
end % function