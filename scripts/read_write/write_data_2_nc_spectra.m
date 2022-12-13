function write_data_2_nc_spectra(data,outfile, config)
% function to write Doppler spectra and related metadata into netcdf4 file
% RG 1.6.2022


%% adjust spectra variables for move convenient use of data

range_offs = data.range_offsets;
range_offs = [range_offs data.n_levels];


for ch = data.no_chirp_seq:-1:1

    % call function to get indexing and new doppler velocity array to shift
    % the Doppler spectra in uniform velocity grid
    [start_ind.(['chirp_' num2str(ch)]), velocity.(['chirp_' num2str(ch)]) ] = ...
        shift_spectra(data.spec(:,range_offs(ch):range_offs(ch+1)-1,1:data.DoppLen(ch)), ...
                      data.velocity(ch,1:data.DoppLen(ch)), ...
                      data.MinVel(:,range_offs(ch):range_offs(ch+1)-1), ...
                      data.MinVel_Correction(:,range_offs(ch):range_offs(ch+1)-1), ...
                      data.DoppMax(ch));


end

%% ################## Create a netCDF file.

ncid = netcdf.create(outfile,'NETCDF4'); 


% ################# Define dimensions
did_time   = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
did_height  = netcdf.defDim(ncid,'height',data.n_levels);
did_no_seq = netcdf.defDim(ncid,'chirp_sequence',data.no_chirp_seq);

% dimensions needed for each chirp
did_vel = zeros(1,data.no_chirp_seq);
for ch = 1:data.no_chirp_seq

    did_vel(ch)    = netcdf.defDim(ncid,['vel_ch' num2str(ch)],length(velocity.(['chirp_' num2str(ch)])));
end


did_hght = zeros(1,data.no_chirp_seq);
n_levels_ch = zeros(1,data.no_chirp_seq);
for ch = 1:data.no_chirp_seq
    n_levels_ch(ch) = range_offs(ch+1)-range_offs(ch);
    did_hght(ch)    = netcdf.defDim(ncid,['height_ch' num2str(ch)], n_levels_ch(ch) );
end

% did_vel(ch)    = netcdf.defDim(ncid,['vel_ch' num2str(ch)],data.DoppLen(ch));
% did_vel    = netcdf.defDim(ncid,'velocity',max(data.DoppLen));

%% ######################## add global attributes

write_globatt_ac3(ncid, config, data, 2)


%% ################ get variable ids and add attributes

% collection of functions defining variables with attributes
defh = outvarmeta;

%%%%%%%%%% coordinate variables

id_time = defh.time(ncid, did_time, 'NC_DOUBLE', isfield(data, 'totsampchangelabel' ));

id_height = defh.height(ncid, did_height);

% not coordinate variables, but following the same order as in moments file
id_lat = defh.lat(ncid);
id_lon = defh.lon(ncid);

% % % 

id_no_seq = defh.chirp_sequence(ncid, did_no_seq);

id_vel = zeros(1,data.no_chirp_seq);
for ch = 1:data.no_chirp_seq
    id_vel(ch) = netcdf.defVar(ncid,['vel_ch' num2str(ch)],'nc_float',did_vel(ch));
    netcdf.putAtt(ncid,id_vel(ch),'long_name',['Doppler velocity array for chirp ' num2str(ch)]);
    netcdf.putAtt(ncid,id_vel(ch),'units','m s-1');
    netcdf.putAtt(ncid,id_vel(ch),'comment',...
        ['The velocity array is asymmetric, i.e. the absolute values of maximum and minumum velocities are not equal. ',...
        'Since the spectrum at -v_nyquist and +v_nyquist is the same, the entry at +v_nyquist was cut. '...
        'If data has been dealiased, it is possible that velocity arrays in different files differ.']);
end

id_hght = zeros(1,data.no_chirp_seq);
for ch = 1:data.no_chirp_seq
    id_hght(ch) = netcdf.defVar(ncid,['height_ch' num2str(ch)],'nc_float',did_hght(ch));
    netcdf.putAtt(ncid,id_hght(ch),'long_name',['height for chirp ' num2str(ch)]);
    netcdf.putAtt(ncid,id_hght(ch),'units','m');
    netcdf.putAtt(ncid,id_hght(ch),'positive','up'); % recommended by CF convention
    netcdf.putAtt(ncid,id_hght(ch),'comment','calculated as range + instrument_altitude');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% variables for spectra %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

id_spec = zeros(1,data.no_chirp_seq);

for ch = 1:data.no_chirp_seq
    
    id_spec(ch) = netcdf.defVar(ncid,['spec_ch' num2str(ch)],'nc_float',[did_vel(ch),did_hght(ch),did_time]);
    netcdf.putAtt(ncid,id_spec(ch),'long_name',['Doppler spectra for chirp ' num2str(ch)]);
    netcdf.putAtt(ncid,id_spec(ch),'units','dB');  
    netcdf.putAtt(ncid,id_spec(ch),'ancillary_variables','quality_flag, ze_calibration');
    netcdf.defVarFill(ncid,id_spec(ch),false,NaN('single'))
    defh.ze_comment(data, ncid, id_spec(ch)) % add comment about ze corrections
end



%%%%%%%% chirp_seq dependent variables
id_range_offsets = defh.range_offsets(ncid,did_no_seq);

id_DoppMax = defh.DoppMax(ncid, did_no_seq);

id_DoppLen = netcdf.defVar(ncid,'spec_length','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_DoppLen,'long_name','number of bins in Doppler spectra of each chirp sequence');
netcdf.defVarFill(ncid,id_DoppLen,false,int32(-999))
netcdf.putAtt(ncid,id_DoppLen,'comment','same as the Doppler FFT; the dimension of the corresponding velocity array might differ from spec_length due to the dealiasing algorithm');


id_nAvg = netcdf.defVar(ncid,'num_avg_spec','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_nAvg,'long_name','number of spectra averaged for each chirp');
netcdf.defVarFill(ncid,id_nAvg,false,int32(-999))
netcdf.putAtt(ncid,id_nAvg,'comment','calculated as Number of averaged chirps within a chirp sequence/spec_length')


id_SeqIntTime = defh.SeqIntTime(ncid, did_no_seq);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% multi-D variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

id_VNoisePow_mean = netcdf.defVar(ncid,'mean_noise','nc_float',[did_height,did_time]);
netcdf.putAtt(ncid,id_VNoisePow_mean,'long_name','Doppler spectrum mean noise');
netcdf.putAtt(ncid,id_VNoisePow_mean,'units','dB');
netcdf.defVarFill(ncid,id_VNoisePow_mean,false,NaN('single'))
netcdf.putAtt(ncid,id_VNoisePow_mean,'comment','calculated from the Doppler spectra following Hildebrand and Sekhon, 1974');

id_VNoisePow_peak = netcdf.defVar(ncid,'peak_noise_pow','nc_float',[did_height,did_time]);
netcdf.putAtt(ncid,id_VNoisePow_peak,'long_name','Doppler spectrum peak noise');
netcdf.putAtt(ncid,id_VNoisePow_peak,'units','dB');
netcdf.defVarFill(ncid,id_VNoisePow_peak,false,NaN('single'))
netcdf.putAtt(ncid,id_VNoisePow_peak,'comment','calculated from the Doppler spectra following Hildebrand and Sekhon, 1974');

if isfield(data, 'SLv')
    id_SLv = defh.SLv(ncid, did_height, did_time);
end
if isfield(data, 'std_noise') % from RPG software version 1
    id_NStd = defh.noisestd(ncid, did_no_seq, did_time);
end

%%%%%%%%%% other variables

id_QF = defh.aggregFlag(ncid, did_time); % quality flags - using same as in moment file

% create here as empty variable for later use
id_ZeCalib = defh.zecalib(ncid);

id_swv = defh.radar_software(ncid);


%% ###################### compression

% % chirp_seq dependent variables
netcdf.defVarDeflate(ncid,id_DoppMax,true,true,9);
netcdf.defVarDeflate(ncid,id_DoppLen,true,true,9);
netcdf.defVarDeflate(ncid,id_nAvg,true,true,9);
netcdf.defVarDeflate(ncid,id_SeqIntTime,true,true,9);
netcdf.defVarDeflate(ncid,id_range_offsets,true,true,9);

% time dependend variables
netcdf.defVarDeflate(ncid,id_QF,true,true,9);

% multi-D variables
netcdf.defVarDeflate(ncid,id_VNoisePow_mean,true,true,9);
netcdf.defVarDeflate(ncid,id_VNoisePow_peak,true,true,9);

if isfield(data, 'SLv')
    netcdf.defVarDeflate(ncid,id_SLv,true,true,9);
end
if isfield(data, 'std_noise') % from RPG software version 1
    netcdf.defVarDeflate(ncid,id_NStd,true,true,9);
end 

% variables for spectra
for ch = 1:data.no_chirp_seq
    netcdf.defVarDeflate(ncid,id_spec(ch),true,true,9);
end

if data.DualPol > 0
    disp('WARMING!!! No polarimetric variables included in the output files')
end


netcdf.endDef(ncid);

%% ####################### put variables into file

% variables for dimensions
netcdf.putVar(ncid,id_no_seq,0,data.no_chirp_seq,1:data.no_chirp_seq);

% scalar variables
netcdf.putVar(ncid,id_lat,data.Lat);
netcdf.putVar(ncid,id_lon,data.Lon);
netcdf.putVar(ncid,id_swv,str2double(data.radarsw));

% range dependet
% netcdf.putVar(ncid,id_range,0,data.n_levels,data.range); DON'T INCLUDE
netcdf.putVar(ncid,id_height,0,data.n_levels,data.range+config.MSL);

% chrip seq dependent variables
netcdf.putVar(ncid,id_DoppMax,0,data.no_chirp_seq,data.DoppMax);
netcdf.putVar(ncid,id_DoppLen,0,data.no_chirp_seq,data.DoppLen);
netcdf.putVar(ncid,id_nAvg,0,data.no_chirp_seq,data.nAvg);
netcdf.putVar(ncid,id_SeqIntTime,0,data.no_chirp_seq,data.SeqIntTime);
netcdf.putVar(ncid,id_range_offsets,0,data.no_chirp_seq,data.range_offsets-1);

% time dependent variables
netcdf.putVar(ncid,id_time,0,data.totsamp, double(data.time) + double(data.sampleTms).*1e-3);

flag_aggregate = ac3_aggregate_flag(data);
netcdf.putVar(ncid,id_QF,0,data.totsamp,flag_aggregate);

% % multidimensional variables

netcdf.putVar(ncid,id_VNoisePow_mean,[0,0],[data.n_levels,data.totsamp],10.*log10(data.VNoisePow_mean'));
netcdf.putVar(ncid,id_VNoisePow_peak,[0,0],[data.n_levels,data.totsamp],10.*log10(data.VNoisePow_peak'));
if isfield(data, 'SLv')
    netcdf.putVar(ncid,id_SLv,[0,0],[data.n_levels,data.totsamp],10.*log10(data.SLv'));
end
if isfield(data, 'std_noise') % from RPG software version 1
    netcdf.putVar(ncid,id_NStd,[0,0],[data.no_chirp_seq,data.totsamp],(data.std_noise'));
end

% variables for spectra


for ch = 1:data.no_chirp_seq
    netcdf.putVar(ncid,id_vel(ch),0,length(velocity.(['chirp_' num2str(ch)])),velocity.(['chirp_' num2str(ch)]));
    netcdf.putVar(ncid,id_hght(ch),0,range_offs(ch+1)-range_offs(ch),data.range(range_offs(ch):range_offs(ch+1)-1)+config.MSL);
end



for ch = data.no_chirp_seq:-1:1
  
    % create new spectra array with new velocity dimensions
    spec_temp = NaN([length(data.time),  range_offs(ch+1) - range_offs(ch), length(velocity.(['chirp_' num2str(ch)]))], 'single');

    % sample the shifted spectra into a fixed velocity dimension using linear indexing
    [sx, sy, sz] = size(spec_temp);
      
    I1 = repmat([1:sx], [1 sy*data.DoppLen(ch)])';
    I2 = reshape(repmat(1:sy, [sx data.DoppLen(ch)]), 1, [])';
    I3 = repmat(start_ind.(['chirp_' num2str(ch)])(:), 1, single(data.DoppLen(ch))) + repmat([0:single(data.DoppLen(ch))-1], sx*sy,1);

    idx = sub2ind([sx, sy, sz], I1, I2, I3(:));
    
    % add shifted spectra in new array
    spec_temp(idx) = data.spec(:,range_offs(ch):range_offs(ch+1)-1,1:data.DoppLen(ch));
    
    netcdf.putVar(ncid,id_spec(ch),[0,0,0],[length(velocity.(['chirp_' num2str(ch)])),n_levels_ch(ch),data.totsamp],permute( 10.*log10(spec_temp) ,[3,2,1]));

    clear spec_temp
end

% netcdf.putVar(ncid,id_vel,[0,0],[max(data.DoppLen),data.no_chirp_seq],data.velocity');
% netcdf.putVar(ncid,id_spec,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute( 10.*log10(data.spec) ,[3,2,1]));


netcdf.close(ncid);


end % function


function [start_ind, vel_out] = shift_spectra(this_spectra, this_vel, this_minvel, this_minvelcorr, this_vnyq)

    % Doppler resolution
    dv = this_vel(2) - this_vel(1);

    % check that minvel exist everywhere were spectra exists
    ind_data = any(~isnan(this_spectra),3);
    
    % replace nan's in minvel, when data for the bin exists - NOTE! this is
    % a workaround, minvel should not be nan and the underlying problem in
    % the dealiasing algorithm should be fixed
    ind_nan = isnan(this_minvel) & ind_data;
    this_minvel(ind_nan) = this_vel(1);
    this_minvelcorr(ind_nan) = 0;
    
    % remove minvel_correction from minvel to get indexing
    minvel_temp = this_minvel - this_minvelcorr;

    ind_shift = (this_vel(1) - minvel_temp)./dv; % number of bins the spectra is shifted, positive shift means minvel < this_vel(1)
    
    % add shift from minvel_correction
    ind_shift = ind_shift - this_minvelcorr/this_vnyq/2 * length(this_vel);
  
    % creating new velocity array
    vel_out = [ fliplr( this_vel(1):-dv:this_vel(1)-max(ind_shift(:))*dv ) this_vel(2:end) ];
    vel_out = [vel_out(1:end-1) vel_out(end):dv:(vel_out(end)-min(ind_shift(:))*dv) ];
    
    % get index in new velocity array
    start_ind = max(ind_shift(:)) - ind_shift+1;
    
    start_ind( isnan(start_ind) ) = 1; % for empty bins, but indexes need to be real numbers
    
end % function
