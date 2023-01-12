function fh = write_data_2_nc_moments(data, outfile, config)
% function to write moments and most important metadata into netcdf4 file
% RG 1.6.2022

%################## Create a netCDF file.

ncid = netcdf.create(outfile,'NETCDF4'); 

%% ################# Define dimensions
did_time   = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
did_height  = netcdf.defDim(ncid,'height',data.n_levels);
did_no_seq = netcdf.defDim(ncid,'chirp_sequence',data.no_chirp_seq);

%% ######################## add global attributes

write_globatt_ac3(ncid, config, data, 1)

if data.cal_mom == 3 % if moments not calculated by our script, add note
    netcdf.putAtt(ncid,glob,'note_dataprocessing', 'Moments are taken from RPG processing sofware (lv1)');
end


%% ################ get variable ids and add attributes

% collection of functions defining variables with attributes
defh = outvarmeta;

%%%%%%%%%% coordinate variables

% id_range = netcdf.defVar(ncid,'range','nc_float',did_height);
% netcdf.putAtt(ncid,id_range,'long_name',['range from the radar antenna to the '...
%                                         'center of the radar range gate']);
% netcdf.putAtt(ncid,id_range,'units','m');
% netcdf.putAtt(ncid,id_range,'positive','up'); % recommended by CF convention


id_time = defh.time(ncid, did_time, 'NC_DOUBLE', isfield(data, 'totsampchangelabel' ));

id_height = defh.height(ncid, did_height);

id_no_seq = defh.chirp_sequence(ncid, did_no_seq);

id_lat = defh.lat(ncid);
id_lon = defh.lon(ncid);

id_MSL = netcdf.defVar(ncid,'instrument_altitude','nc_float',[]);
netcdf.putAtt(ncid,id_MSL,'long_name','instrument altitude above mean sea level');
netcdf.putAtt(ncid,id_MSL,'units','m');
netcdf.defVarFill(ncid,id_MSL,false,NaN('single'))

%%%%%%%%%% other variables in a convenient order for reading file headers

%--- radar moments ---

id_Ze = netcdf.defVar(ncid,'ze','nc_float',[did_height,did_time]);
netcdf.putAtt(ncid,id_Ze,'long_name','equivalent radar reflectivity factor');
netcdf.putAtt(ncid,id_Ze,'standard_name','equivalent_reflectivity_factor');
netcdf.putAtt(ncid,id_Ze,'units','dB');
netcdf.putAtt(ncid,id_Ze,'ancillary_variables','quality_flag, ze_calibration');
netcdf.defVarFill(ncid,id_Ze,false,NaN('single'))
defh.ze_comment(data, ncid, id_Ze) % add comment about ze corrections

id_vm = netcdf.defVar(ncid,'vm','nc_float',[did_height,did_time]);
netcdf.putAtt(ncid,id_vm,'long_name','mean Doppler velocity');
netcdf.putAtt(ncid,id_vm,'units','m s-1');
netcdf.putAtt(ncid,id_vm,'ancillary_variables','quality_flag');    
netcdf.defVarFill(ncid,id_vm,false,NaN('single'))
netcdf.putAtt(ncid,id_vm,'comment',['negative velocities indicate particles moving downwards'])

id_sigma = netcdf.defVar(ncid,'sw','nc_float',[did_height,did_time]);
netcdf.putAtt(ncid,id_sigma,'long_name','Doppler spectrum width');
netcdf.putAtt(ncid,id_sigma,'units','m s-1');
netcdf.putAtt(ncid,id_sigma,'ancillary_variables','quality_flag');    
netcdf.defVarFill(ncid,id_sigma,false,NaN('single'))

id_skew = netcdf.defVar(ncid,'skew','nc_float',[did_height,did_time]);
netcdf.putAtt(ncid,id_skew,'long_name','Doppler spectrum skewness');
netcdf.putAtt(ncid,id_skew,'units','unitless');
netcdf.putAtt(ncid,id_skew,'ancillary_variables','quality_flag');
netcdf.defVarFill(ncid,id_skew,false,NaN('single'))

id_kurt = netcdf.defVar(ncid,'kurt','nc_float',[did_height,did_time]);
netcdf.putAtt(ncid,id_kurt,'long_name','Doppler spectrum kurtosis');
netcdf.putAtt(ncid,id_kurt,'units','unitless');
netcdf.putAtt(ncid,id_kurt,'ancillary_variables','quality_flag');
netcdf.defVarFill(ncid,id_kurt,false,NaN('single'))

if data.DualPol > 0
    disp('WARMING!!! No polarimetric variables included in the output files')
end

if isfield(data, 'SLv')
    id_SLv = defh.SLv(ncid, did_height, did_time);
end

if isfield(data, 'std_noise') % from RPG software version 1
    id_NStd = defh.noisestd(ncid, did_no_seq, did_time);
end

if isfield(data, 'QF')
    id_QF = defh.aggregFlag(ncid, did_time, true);
else
    id_QF = defh.aggregFlag(ncid, did_time, false);
end


% create here as empty variable for later use
id_ZeCalib = defh.zecalib(ncid);

% technical parameters for "expert use"
id_freq = netcdf.defVar(ncid,'frequency','nc_float',[]);
netcdf.putAtt(ncid,id_freq,'long_name','central transmission frequency');
netcdf.putAtt(ncid,id_freq,'standard_name','radiation_frequency');
netcdf.putAtt(ncid,id_freq,'units','GHz');
netcdf.defVarFill(ncid,id_freq,false,NaN('single'))

id_HPBW = netcdf.defVar(ncid,'beam_width','nc_float',[]);
netcdf.putAtt(ncid,id_HPBW,'long_name','antenna half power beam width');
netcdf.putAtt(ncid,id_HPBW,'units','degrees');
netcdf.defVarFill(ncid,id_HPBW,false,NaN('single'))

id_CalInt = netcdf.defVar(ncid,'zerocal_interval','nc_int',[]);
netcdf.putAtt(ncid,id_CalInt,'long_name','sample interval (number of samples) between automated zero calibrations');
netcdf.putAtt(ncid,id_CalInt,'units','count');
netcdf.defVarFill(ncid,id_CalInt,false,int32(-999))

id_AntG = netcdf.defVar(ncid,'antenna_gain','nc_float',[]);
netcdf.putAtt(ncid,id_AntG,'long_name','linear antenna gain');
netcdf.putAtt(ncid,id_AntG,'units','unitless');
netcdf.defVarFill(ncid,id_AntG,false,NaN('single'))

id_swv = defh.radar_software(ncid);

id_DoppMax = defh.DoppMax(ncid, did_no_seq);

id_range_offsets = defh.range_offsets(ncid,did_no_seq);
                                           
id_SeqIntTime = defh.SeqIntTime(ncid, did_no_seq);


id_blowstatus = netcdf.defVar(ncid,'blower_status','nc_byte',did_time);
netcdf.putAtt(ncid,id_blowstatus,'long_name','blower_status');
netcdf.putAtt(ncid,id_blowstatus,'standard_name','status_flag');
netcdf.putAtt(ncid,id_blowstatus,'flag_values', [0b0, 0b1]);
netcdf.putAtt(ncid,id_blowstatus,'flag_meanings','blower_off blower_on');
netcdf.defVarFill(ncid,id_blowstatus,false,0b10)


id_heatstatus = netcdf.defVar(ncid,'heater_status','nc_byte',did_time);
netcdf.putAtt(ncid,id_heatstatus,'long_name','heater_status');
netcdf.putAtt(ncid,id_heatstatus,'standard_name','status_flag');
netcdf.putAtt(ncid,id_heatstatus,'flag_values', [0b0, 0b1]);
netcdf.putAtt(ncid,id_heatstatus,'flag_meanings','heater_off heater_on');
netcdf.defVarFill(ncid,id_heatstatus,false,0b10)


id_TransPow = netcdf.defVar(ncid,'transmitter_power','nc_float',did_time);
netcdf.putAtt(ncid,id_TransPow,'long_name','transmitter_power');
netcdf.putAtt(ncid,id_TransPow,'units','W');
netcdf.defVarFill(ncid,id_TransPow,false,NaN('single'))





%% ###################### compression

% chirp_seq dependent variables
netcdf.defVarDeflate(ncid,id_DoppMax,true,true,9);
netcdf.defVarDeflate(ncid,id_range_offsets,true,true,9);
netcdf.defVarDeflate(ncid,id_SeqIntTime,true,true,9);

% time dependent variables
netcdf.defVarDeflate(ncid,id_blowstatus,true,true,9);
netcdf.defVarDeflate(ncid,id_heatstatus,true,true,9);
netcdf.defVarDeflate(ncid,id_TransPow,true,true,9);
netcdf.defVarDeflate(ncid,id_QF,true,true,9);

% multi-D variables
netcdf.defVarDeflate(ncid,id_Ze,true,true,9);
netcdf.defVarDeflate(ncid,id_vm,true,true,9);
netcdf.defVarDeflate(ncid,id_sigma,true,true,9);
netcdf.defVarDeflate(ncid,id_skew,true,true,9);
netcdf.defVarDeflate(ncid,id_kurt,true,true,9);

if isfield(data, 'SLv')
    netcdf.defVarDeflate(ncid,id_SLv,true,true,9);
end
if isfield(data, 'std_noise') % from RPG software version 1
    netcdf.defVarDeflate(ncid,id_NStd,true,true,9);
end 

netcdf.endDef(ncid);

%% ####################### put variables into file

% variables for dimensions
netcdf.putVar(ncid,id_no_seq,0,data.no_chirp_seq,1:data.no_chirp_seq);

% scalars
netcdf.putVar(ncid,id_lat,data.Lat);
netcdf.putVar(ncid,id_lon,data.Lon);
netcdf.putVar(ncid,id_MSL,config.MSL);
netcdf.putVar(ncid,id_freq,data.freq);
netcdf.putVar(ncid,id_HPBW,data.HPBW);
netcdf.putVar(ncid,id_CalInt,data.CalInt);
netcdf.putVar(ncid,id_AntG,data.AntG);
netcdf.putVar(ncid,id_swv,str2num(data.radarsw));


% range dependent
netcdf.putVar(ncid,id_height,0,data.n_levels,data.range+config.MSL);


% chrip seq dependent variables
netcdf.putVar(ncid,id_DoppMax,0,data.no_chirp_seq,data.DoppMax);
netcdf.putVar(ncid,id_range_offsets,0,data.no_chirp_seq,data.range_offsets-1);
netcdf.putVar(ncid,id_SeqIntTime,0,data.no_chirp_seq,data.SeqIntTime);

% time dependent variables
netcdf.putVar(ncid,id_time,0,data.totsamp,  double(data.time) + double(data.sampleTms).*1e-3);
netcdf.putVar(ncid,id_TransPow,0,data.totsamp,data.TransPow);

% convert to have 2 for missing value
data.blower(isnan(data.blower)) = 2;
netcdf.putVar(ncid,id_blowstatus,0,data.totsamp,data.blower);

% convert to have 2 for missing value
data.heater(isnan(data.heater)) = 2;
netcdf.putVar(ncid,id_heatstatus,0,data.totsamp,data.heater);

flag_aggregate = ac3_aggregate_flag(data);
netcdf.putVar(ncid,id_QF,0,data.totsamp,flag_aggregate);

% multidimensional variables
netcdf.putVar(ncid,id_Ze,[0,0],[data.n_levels,data.totsamp],10.*log10(data.Ze'));
netcdf.putVar(ncid,id_vm,[0,0],[data.n_levels,data.totsamp],data.vm');
netcdf.putVar(ncid,id_sigma,[0,0],[data.n_levels,data.totsamp],data.sigma');
netcdf.putVar(ncid,id_skew,[0,0],[data.n_levels,data.totsamp],data.skew');
netcdf.putVar(ncid,id_kurt,[0,0],[data.n_levels,data.totsamp],data.kurt');

if isfield(data, 'SLv')
    netcdf.putVar(ncid,id_SLv,[0,0],[data.n_levels,data.totsamp],data.SLv');
end

if isfield(data, 'std_noise') % from RPG software version 1
    netcdf.putVar(ncid,id_NStd,[0,0],[data.no_chirp_seq,data.totsamp],(data.std_noise'));
end

netcdf.close(ncid);

end % function

