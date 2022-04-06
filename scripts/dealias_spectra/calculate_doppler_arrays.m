function velocitymatrix = calculate_doppler_arrays(data)
% function to calculate the Doppler velocity array for each time-range bin,
% following the dealiasing algorithm that might have shifted these
% Rosa Gierens 7.6.2021

velocitymatrix = NaN(size(data.spec));

for i = 1:data.totsamp
    for j = 1:data.n_levels
        
        % find number of current chirp index
        chirp_idx = int32(find(data.range_offsets(2:end) - j > 0,1,'first'));
        if isempty(chirp_idx) % then j is within last chirp
            chirp_idx = data.no_chirp_seq;
        end
    
        if isnan( data.MinVel(i,j)) % not sure if this is necessary
            velocitymatrix(i,j,:) = data.velocity(chirp_idx,:);
            
        else
            % recontsruct velocity array for each bin
            % vel_true = velocity + Minvel(i,j) - velocity(1);

            velocitymatrix(i,j,:) = data.velocity(chirp_idx,:) +  data.MinVel(i,j) - data.velocity(chirp_idx,1);
        end          
    end 
end
