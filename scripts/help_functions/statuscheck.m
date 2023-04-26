function data = statuscheck(data)
% found some cases where status got unreasonable values, these happened in 
% files where overwrote time with lv1 file timestamp, so problaly bad data
% logged and not cathed by time check as before.
% This function checks for timestamps where status gets any value that is
% not 0, 1, 10 or 11, and removes all data related to that time stamp from 
% the dataset
% RG 11.4.2023



% check if in any profile reflectivity > 50 dB 
if all(data.status == 0 | data.status == 1 | data.status == 10 | data.status == 11)
    return
    
else

    % find profile with bad data
    indout = ~(data.status == 0 | data.status == 1 | data.status == 10 | data.status == 11);


    % loop over fieldnames of data
    fn = fieldnames(data);
    for k=1:numel(fn)

        if any(size(data.(fn{k})) == data.totsamp) % check if time dimenstion

            timdim = find(size(data.(fn{k})) == data.totsamp);


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
    
    %data.totsamp = data.totsamp-1; % 
end

% making sure that not totsamp is the right size
data.totsamp = length(data.time);