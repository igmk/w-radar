function data = zesanitycheck(data)
% Sometimes Ze gets unphysical values, this seems to happen sometimes 
% when the measurement is interrupted in the last logged profile
% This function checks for profiles where spectral reflectivity > 50 dB,
% and removes entire profile from the dataset
% RG 7.9.2021


% check if in any profile reflectivity > 50 dB 
if ~any(data.spec(:) > 10^(50/10))
    return
    
else

    % find profile with bad data
    indout = any(any(data.spec > 10^(50/10), 3)');


    % loop over fieldnames of data
    fn = fieldnames(data);
    for k=1:numel(fn)

        if any(size(data.(fn{k})) == data.totsamp) % check if time dimenstion

            timdim = find(size(data.(fn{k})) == data.totsamp);
            
            % in this file, at this point range and time dimensions have same length!
            if data.starttime(1) == 593434800 
                if length(timdim) > 1
                    timdim = timdim(1); 
                end
                
                if strcmp(fn{k}, 'range') || strcmp(fn{k}, 'Fr')
                    continue
                end
            end
            

            switch length(size(data.(fn{k}))) % remove bad entry on time dimensions

                case 1
                    data.(fn{k})(indout) = [];

                case 2
                     switch timdim
                         case 1
                             data.(fn{k})(indout,:) = [];

                         case 2
                             data.(fn{k})(:,indout) = [];
                     end

                case 3

                    switch timdim
                        case 1
                            data.(fn{k})(indout,:,:) = [];
                        case 2
                            data.(fn{k})(:,indout,:) = [];
                        case 3
                            data.(fn{k})(:,:,indout) = [];
                    end
            end


        end

    end
    
    data.totsamp = data.totsamp-1; % 
end
