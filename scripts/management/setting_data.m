function [data] = setting_data(data, config, path)

specsize = size(data.spec);


% ############## add some variables
data_fieldnames = fieldnames(data);

% number of spectral averages
data.nAvg = data.SeqAvg./data.DoppLen; % number of spectral averages

% create velocity array
dv = 2.*data.DoppMax./double(data.DoppLen);
data.velocity = NaN(data.no_chirp_seq,max(data.DoppLen));
for i = 1:data.no_chirp_seq
    data.velocity(i,1:data.DoppLen(i)) = -data.DoppMax(i):dv(i):data.DoppMax(i)-dv(i);
end

%minimum velocity 
if ~any(strcmp(data_fieldnames,'MinVel'))
    data.MinVel = zeros(specsize(1),specsize(2)); % contains minimum velocity in spectra after dealiasing was applied
else
    data.MinVel(data.MinVel == -999) = 0;
end

%Dealising applied?
if ~any(strcmp(data_fieldnames,'Aliasmask'))
    data.Aliasmask = NaN(specsize(1),specsize(2)); %indicating if dealiasing was applied
end

% NoisePow not provided? given if CompEna > 0
if ~any(strcmp(data_fieldnames,'VNoisePow_mean')) % then NoisePow is not provided by RPG
    data.VNoisePow_mean = NaN(specsize(1),specsize(2)); % preallocate mean noise floor
else
    data.VNoisePow_mean(data.VNoisePow_mean == -999) = NaN;
    for i = 1:data.no_chirp_seq % convert into bin noise power            
        if i == data.no_chirp_seq
            r_idx = data.range_offsets(i):specsize(2);
        else
            r_idx = data.range_offsets(i):data.range_offsets(i+1);
        end            
        data.VNoisePow_mean(:,r_idx) = data.VNoisePow_mean(:,r_idx)./single(data.DoppLen(i));
        if data.DualPol > 0                
            data.HNoisePow_mean(data.HNoisePow_mean == -999) = NaN;
            data.HNoisePow_mean(:,r_idx) = data.HNoisePow_mean(:,r_idx)./single(data.DoppLen(i));
        end            
    end
end

% preallocate peak noise floor
data.VNoisePow_peak = NaN(specsize(1),specsize(2)); 

% preallocate peak noise floor
if data.DualPol > 0
    data.HNoisePow_peak = NaN(specsize(1),specsize(2)); 
end

% create dummy status flag
data.AliasStatus = NaN(specsize(1),specsize(2)); % dummy variable that provides information on quality of dealiasing, not given by RPG software

%Checking if dealias was already applied by RPG software.
if ~data.AntiAlias  
    
    % set MinVel matrix
    
    % if dealiasing was performed then the true velocity array can be
    % calculated as the following:
    % vel_true =  data.velocity + data.Minvel(i,j) - data.velocity(1)
    % note that the velocity array is different in every chrip sequence
    for ii = 1:data.no_chirp_seq % set MinVel to first entry of velocity
        if ii == data.no_chirp_seq
            r_idx = data.range_offsets(ii):specsize(2);
        else
            r_idx = data.range_offsets(ii):data.range_offsets(ii+1);
        end
        data.MinVel(:,r_idx) = data.velocity(ii,1);
    end
end

% create correction matrix for data.MinVel
% in dealias_spectra_vm_column_quality_check.m wrongly dealiased
% velocity bins (i.e. shifted by k*2*v_n) are corrected. the corrected
% offset will be stored in data.MinVel_Correction. data.MinVel is
% stored offset corrected
% will only be modified if dealiasing is performed
data.MinVel_Correction = zeros(specsize(1),specsize(2));    


% flag for compressed spectra
if (isfield(data, 'CompEna') && ne(data.CompEna,0))
    data.compress_spec = true;
else
    data.compress_spec = false;
end

% create variable for chirp integration times
data.ChirpIntTime = NaN(1,data.no_chirp_seq);

% calculate chirp integration times
data.ChirpIntTime = ChirpIntegrationTimes(data.DoppMax, data.freq, data.SeqAvg);


% set -999 values to NaN
data = fill2nan_struct(data,-999); 
    

% separate blower and heater flags 
data.blower = NaN(size(data.status));
data.heater = NaN(size(data.status));

% 0/1 = heater on/off; 0/10 = blower on/off.

% set blower status
data.blower( data.status >= 10 ) = 1;
data.blower( data.status < 10 ) = 0;

% heater status is what is left when removing blower status
data.heater = data.status;
data.heater(data.status >= 10) = data.heater(data.status >= 10) - 10;

% store for later
data.infullpath = path.lv0; 


