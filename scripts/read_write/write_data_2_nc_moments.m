function fh = write_data_2_nc_moments(data, outfile, config)
% function to write moments and most important metadata into netcdf4 file
% RG 1.6.2022

%################## Create a netCDF file.

ncid = netcdf.create(outfile,'NETCDF4'); 

%% ################# Define dimensions
did_time   = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
did_height  = netcdf.defDim(ncid,'height',data.n_levels);

%% ######################## add global attributes

write_globatt_ac3(ncid, config, data, 1)

if data.cal_mom == 3 % if moments not calculated by our script, add note
    netcdf.putAtt(ncid,glob,'note_dataprocessing', 'Moments are taken from RPG processing sofware (lv1)');
end

%%%%%%%%%% create group for technical variables

ncid_tech = netcdf.defGrp(ncid,'Technical_parameters');

if isfield(data, 'std_noise') % this variable included in the global group, from RPG software version 1
    did_no_seq = netcdf.defDim(ncid,'chirp_sequence',data.no_chirp_seq);
else % standard set-up
    did_no_seq = netcdf.defDim(ncid_tech,'chirp_sequence',data.no_chirp_seq);
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


id_time = defh.time(ncid, did_time, isfield(data, 'totsampchangelabel' ));

id_height = defh.height(ncid, did_height);

id_lat = defh.lat(ncid);
id_lon = defh.lon(ncid);


%%%%%%%%%% other variables in a convenient order for reading file headers
id_sampleTms = defh.sampleTms(ncid, did_time);

%--- radar moments ---

id_Ze = netcdf.defVar(ncid,'ze','nc_float',[did_height,did_time]);
netcdf.putAtt(ncid,id_Ze,'long_name','equivalent radar reflectivity factor');
netcdf.putAtt(ncid,id_Ze,'standard_name','equivalent_reflectivity_factor');
netcdf.putAtt(ncid,id_Ze,'units','dBz');
netcdf.putAtt(ncid,id_Ze,'ancillary_variables','quality_flag');
netcdf.putAtt(ncid,id_Ze,'_FillValue',NaN('single'));
if isfield(data, 'Ze_label') % Ze corrected, adding note
    netcdf.putAtt(ncid,id_Ze,'comment',data.Ze_label);
    netcdf.putAtt(ncid,id_Ze,'corretion_dB',data.Ze_corr);
end

id_vm = netcdf.defVar(ncid,'vm','nc_float',[did_height,did_time]);
netcdf.putAtt(ncid,id_vm,'long_name','mean Doppler velocity');
netcdf.putAtt(ncid,id_vm,'units','m s-1');
netcdf.putAtt(ncid,id_vm,'ancillary_variables','quality_flag');    
netcdf.putAtt(ncid,id_vm,'_FillValue',NaN('single'));
netcdf.putAtt(ncid,id_vm,'comment',['negative velocities indicate particles moving downwards'])

id_sigma = netcdf.defVar(ncid,'sw','nc_float',[did_height,did_time]);
netcdf.putAtt(ncid,id_sigma,'long_name','Doppler spectrum width');
netcdf.putAtt(ncid,id_sigma,'units','m s-1');
netcdf.putAtt(ncid,id_sigma,'ancillary_variables','quality_flag');    
netcdf.putAtt(ncid,id_sigma,'_FillValue',NaN('single'));

id_skew = netcdf.defVar(ncid,'skew','nc_float',[did_height,did_time]);
netcdf.putAtt(ncid,id_skew,'long_name','Doppler spectrum skewness');
netcdf.putAtt(ncid,id_skew,'ancillary_variables','quality_flag');
netcdf.putAtt(ncid,id_skew,'_FillValue',NaN('single'));

id_kurt = netcdf.defVar(ncid,'kurt','nc_float',[did_height,did_time]);
netcdf.putAtt(ncid,id_kurt,'long_name','Doppler spectrum kurtosis');
netcdf.putAtt(ncid,id_kurt,'ancillary_variables','quality_flag');
netcdf.putAtt(ncid,id_kurt,'_FillValue',NaN('single'));

if data.DualPol > 0
    disp('WARMING!!! No polarimetric variables included in the output files')
end

if isfield(data, 'SLv')
    id_SLv = defh.SLv(ncid, did_height, did_time);
end

if isfield(data, 'std_noise') % from RPG software version 1
    id_NStd = defh.noisestd(ncid, did_no_seq, did_time);
end


if isfield(data, 'std_noise') % from RPG software version 1
    id_no_seq = defh.chirp_sequence(ncid, did_no_seq);    
else
    id_no_seq = defh.chirp_sequence(ncid_tech, did_no_seq);
end

% technical parameters for "expert use"

id_MSL = netcdf.defVar(ncid_tech,'instrument_altitude','nc_float',[]);
netcdf.putAtt(ncid_tech,id_MSL,'long_name','instrument altitude above mean sea level');
netcdf.putAtt(ncid_tech,id_MSL,'units','m');
netcdf.putAtt(ncid_tech,id_MSL,'_FillValue',NaN('single'));

id_freq = netcdf.defVar(ncid_tech,'frequency','nc_float',[]);
netcdf.putAtt(ncid_tech,id_freq,'long_name','central transmission frequency');
netcdf.putAtt(ncid_tech,id_freq,'standard_name','radiation_frequency');
netcdf.putAtt(ncid_tech,id_freq,'units','GHz');
netcdf.putAtt(ncid_tech,id_freq,'_FillValue',NaN('single'));

id_HPBW = netcdf.defVar(ncid_tech,'beam_width','nc_float',[]);
netcdf.putAtt(ncid_tech,id_HPBW,'long_name','antenna half power beam width');
netcdf.putAtt(ncid_tech,id_HPBW,'units','degrees');
netcdf.putAtt(ncid_tech,id_HPBW,'_FillValue',NaN('single'));
                                  
id_CalInt = netcdf.defVar(ncid_tech,'zerocal_interval','nc_int',[]);
netcdf.putAtt(ncid_tech,id_CalInt,'long_name','sample interval (number of samples) between automated zero calibrations');
netcdf.putAtt(ncid_tech,id_CalInt,'units','count');
netcdf.putAtt(ncid_tech,id_CalInt,'_FillValue',int32(-999));

id_AntG = netcdf.defVar(ncid_tech,'antenna_gain','nc_float',[]);
netcdf.putAtt(ncid_tech,id_AntG,'long_name','linear antenna gain');
netcdf.putAtt(ncid_tech,id_AntG,'_FillValue',NaN('single'));

id_swv = defh.radar_software(ncid_tech);

id_DoppMax = defh.DoppMax(ncid_tech, did_no_seq);

id_range_offsets = defh.range_offsets(ncid_tech,did_no_seq);
                                           
id_SeqIntTime = defh.SeqIntTime(ncid_tech, did_no_seq);


id_blowstatus = netcdf.defVar(ncid_tech,'blower_status','nc_byte',did_time);
netcdf.putAtt(ncid_tech,id_blowstatus,'long_name','blower_status');
netcdf.putAtt(ncid_tech,id_blowstatus,'flag_values', [0b0, 0b1, 0b10]);
netcdf.putAtt(ncid_tech,id_blowstatus,'flag_meanings','missing_data blower_off blower_on');
netcdf.putAtt(ncid_tech,id_blowstatus,'standard_name','status_flag');


id_heatstatus = netcdf.defVar(ncid_tech,'heater_status','nc_byte',did_time);
netcdf.putAtt(ncid_tech,id_heatstatus,'long_name','heater_status');
netcdf.putAtt(ncid_tech,id_heatstatus,'flag_values', [0b0, 0b1, 0b10]);
netcdf.putAtt(ncid_tech,id_heatstatus,'flag_meanings','missing_data heater_off heater_on');
netcdf.putAtt(ncid_tech,id_heatstatus,'standard_name','status_flag');



id_TransPow = netcdf.defVar(ncid_tech,'transmitter_power','nc_float',did_time);
netcdf.putAtt(ncid_tech,id_TransPow,'long_name','transmitter_power');
netcdf.putAtt(ncid_tech,id_TransPow,'units','W');
netcdf.putAtt(ncid_tech,id_TransPow,'_FillValue',NaN('single'));


% TODO: better flaggin system - this is a placeholder flag                                       
id_QF = netcdf.defVar(ncid,'quality_flag','nc_byte',did_time);
netcdf.putAtt(ncid,id_QF,'standard_name','aggregate_quality_flag');
netcdf.putAtt(ncid,id_QF,'comment','place holder until quality flag redesigned');

% netcdf.putAtt(ncid,id_QF,'long_name',['Quality flag related to the radar operation '... 
%                                       'given by radar software from RPG']);                              
% netcdf.putAtt(ncid,id_QF,'comment', ['To get the bit entries, one has to'...
%                                      'convert the integer into a 4 bit binary. '...
%                                      'bit4 = ADC saturation, bit3 = spectral '...
%                                      'width too high, bit2 = no transmitter '...
%                                      'power leveling. Note that in the above '...
%                                      'convention holds: bit1 = 2^3, '...
%                                      'bit2 = 2^2, bit3 = 2^1, bit4 = 2^0'])




%% ###################### compression

% chirp_seq dependent variables
netcdf.defVarDeflate(ncid_tech,id_DoppMax,true,true,9);
netcdf.defVarDeflate(ncid_tech,id_range_offsets,true,true,9);
netcdf.defVarDeflate(ncid,id_SeqIntTime,true,true,9);

% time dependent variables
netcdf.defVarDeflate(ncid,id_sampleTms,true,true,9);
netcdf.defVarDeflate(ncid_tech,id_blowstatus,true,true,9);
netcdf.defVarDeflate(ncid_tech,id_heatstatus,true,true,9);
netcdf.defVarDeflate(ncid_tech,id_TransPow,true,true,9);
netcdf.defVarDeflate(ncid,id_QF,true,true,9);

% multi-D variables
% netcdf.defVarDeflate(ncid,id_mask,true,true,9);
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

if isfield(data, 'std_noise') % this variable included in the global group, from RPG software version 1
    netcdf.putVar(ncid,id_no_seq,0,data.no_chirp_seq,1:data.no_chirp_seq);
else
    netcdf.putVar(ncid_tech,id_no_seq,0,data.no_chirp_seq,1:data.no_chirp_seq);
end

% scalars
netcdf.putVar(ncid,id_lat,data.Lat);
netcdf.putVar(ncid,id_lon,data.Lon);
netcdf.putVar(ncid_tech,id_MSL,config.MSL);
netcdf.putVar(ncid_tech,id_freq,data.freq);
netcdf.putVar(ncid_tech,id_HPBW,data.HPBW);
netcdf.putVar(ncid_tech,id_CalInt,data.CalInt);
netcdf.putVar(ncid_tech,id_AntG,data.AntG);
netcdf.putVar(ncid_tech,id_swv,str2num(data.radarsw));


% range dependent
netcdf.putVar(ncid,id_height,0,data.n_levels,data.range+config.MSL);


% chrip seq dependent variables
netcdf.putVar(ncid_tech,id_DoppMax,0,data.no_chirp_seq,data.DoppMax);
netcdf.putVar(ncid_tech,id_range_offsets,0,data.no_chirp_seq,data.range_offsets-1);
netcdf.putVar(ncid_tech,id_SeqIntTime,0,data.no_chirp_seq,data.SeqIntTime);

% time dependent variables
netcdf.putVar(ncid,id_time,0,data.totsamp,data.time);
netcdf.putVar(ncid,id_sampleTms,0,data.totsamp,data.sampleTms);
netcdf.putVar(ncid_tech,id_TransPow,0,data.totsamp,data.TransPow);

% convert to have 0 for missing value
data.blower = data.blower + 1;
data.blower(isnan(data.blower)) = 0;
netcdf.putVar(ncid_tech,id_blowstatus,0,data.totsamp,data.blower);

% convert to have 0 for missing value
data.heater = data.heater + 1; 
data.heater(isnan(data.heater)) = 0;
netcdf.putVar(ncid_tech,id_heatstatus,0,data.totsamp,data.heater);

% netcdf.putVar(ncid,id_QF,0,data.totsamp,data.QF);

% multidimensional variables
% netcdf.putVar(ncid,id_mask,[0,0],[data.n_levels,data.totsamp],data.mask');
netcdf.putVar(ncid,id_Ze,[0,0],[data.n_levels,data.totsamp],10.*log10(data.Ze'));
netcdf.putVar(ncid,id_vm,[0,0],[data.n_levels,data.totsamp],data.vm');
netcdf.putVar(ncid,id_sigma,[0,0],[data.n_levels,data.totsamp],data.sigma');
netcdf.putVar(ncid,id_skew,[0,0],[data.n_levels,data.totsamp],data.skew');
netcdf.putVar(ncid,id_kurt,[0,0],[data.n_levels,data.totsamp],data.kurt');

if isfield(data, 'SLv')
    netcdf.putVar(ncid,id_SLv,[0,0],[data.n_levels,data.totsamp],10.*log10(data.SLv'));
end

if isfield(data, 'std_noise') % from RPG software version 1
    netcdf.putVar(ncid,id_NStd,[0,0],[data.no_chirp_seq,data.totsamp],(data.std_noise'));
end

netcdf.close(ncid);

end % function

