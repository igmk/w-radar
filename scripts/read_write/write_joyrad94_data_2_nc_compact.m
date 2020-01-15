function write_joyrad94_data_2_nc_compact(data, outfile, config)

% this function writes joyrad94 data into netcdf4
% Changes of all the 

%% ################## Create a netCDF file.

ncid = netcdf.create(outfile,'NETCDF4'); 


%% ################# Define dimensions

%did_time = netcdf.defDim(ncid,'time',data.totsamp);
%Changed by J.A. Bravo-Aranda
did_time = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
did_range = netcdf.defDim(ncid,'range',data.n_levels);
did_no_seq = netcdf.defDim(ncid,'chirp_sequences',data.no_chirp_seq);
did_scalar = netcdf.defDim(ncid,'scalar',1);


%% ######################## add global attributes
glob = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,glob,'FillValue','NaN');
netcdf.putAtt(ncid,glob,'program_name',data.progname);
if data.modelno == 0
    model = '94 GHz single pol.';
else
    model = '94 GHz dual pol.';
end

netcdf.putAtt(ncid,glob,'PI_NAME',config.pi_name);
netcdf.putAtt(ncid,glob,'PI_AFFILIATION',config.pi_affiliation);
netcdf.putAtt(ncid,glob,'PI_ADDRESS',config.pi_address);
netcdf.putAtt(ncid,glob,'PI_MAIL',config.pi_mail);
netcdf.putAtt(ncid,glob,'DO_NAME',config.do_name);
netcdf.putAtt(ncid,glob,'DO_AFFILIATION',config.do_affiliation);
netcdf.putAtt(ncid,glob,'DO_ADDRESS',config.do_address);
netcdf.putAtt(ncid,glob,'DO_MAIL',config.do_mail);
netcdf.putAtt(ncid,glob,'DS_NAME',config.ds_name);
netcdf.putAtt(ncid,glob,'DS_AFFILIATION',config.ds_affiliation);
netcdf.putAtt(ncid,glob,'DS_ADDRESS',config.ds_address);
netcdf.putAtt(ncid,glob,'DS_MAIL',config.ds_mail);
netcdf.putAtt(ncid,glob,'DATA_DESCRIPTION',config.data_description);
netcdf.putAtt(ncid,glob,'DATA_DISCIPLINE',config.data_discipline);
netcdf.putAtt(ncid,glob,'DATA_GROUP',config.data_group);
netcdf.putAtt(ncid,glob,'DATA_LOCATION',config.data_location);
netcdf.putAtt(ncid,glob,'DATA_SOURCE',config.data_source);
netcdf.putAtt(ncid,glob,'DATA_PROCESSING',config.processing_script);
netcdf.putAtt(ncid,glob,'FILL_VALUE','NaN');
netcdf.putAtt(ncid,glob,'INSTRUMENT_MODEL',model);
netcdf.putAtt(ncid,glob,'MDF_PROGRAM_USED',data.progname);


%% ################ get variable ids and add attributes

%%%%%%%%%% scalar variables

id_lat = netcdf.defVar(ncid,'lat','nc_float',did_scalar);
netcdf.putAtt(ncid,id_lat,'standard_name','LATITUDE.INSTRUMENT');
netcdf.putAtt(ncid,id_lat,'long_name','LATITUDE of measurement site the instrument is located');
netcdf.putAtt(ncid,id_lat,'units','degrees_north');
netcdf.putAtt(ncid,id_lat,'comment','LATITUDE in degrees north [-90,90]');

id_lon = netcdf.defVar(ncid,'lon','nc_float',did_scalar);
netcdf.putAtt(ncid,id_lon,'standard_name','LONGITUDE.INSTRUMENT');
netcdf.putAtt(ncid,id_lon,'long_name','LONGITUDE of measurement site the instrument is located');
netcdf.putAtt(ncid,id_lon,'units','degrees_east');
netcdf.putAtt(ncid,id_lon,'comment','LONGITUDE in degrees east [-180,180]');

id_MSL = netcdf.defVar(ncid,'zsl','nc_float',did_scalar);
netcdf.putAtt(ncid,id_MSL,'standard_name','ALTITUDE.INSTRUMENT');
netcdf.putAtt(ncid,id_MSL,'long_name','ALTITUDE of the measurement site above mean sea level');
netcdf.putAtt(ncid,id_MSL,'units','m');
netcdf.putAtt(ncid,id_MSL,'comment','Height above mean sea level');

id_freq = netcdf.defVar(ncid,'freq_sb','nc_float',did_scalar);
netcdf.putAtt(ncid,id_freq,'standard_name','FREQUENCY');
netcdf.putAtt(ncid,id_freq,'long_name','Transmission FREQUENCY of the radar system');
netcdf.putAtt(ncid,id_freq,'units','s-1');
netcdf.putAtt(ncid,id_freq,'comment','FREQUENCY can be converted into WAVELENGHT = c/FREQUENCY; where c is the speed of light');

id_AntiAlias = netcdf.defVar(ncid,'anti_alias','nc_byte',did_scalar);
netcdf.putAtt(ncid,id_AntiAlias,'standard_name','processing.quality.flag.anti_alias');
netcdf.putAtt(ncid,id_AntiAlias,'long_name','Quality flag for dealiasing');
netcdf.putAtt(ncid,id_AntiAlias,'comment','The falg index shows: 0 = no dealiasing applied, 1 = dealiasing by RPG, 2 = dealiasing by the applied code (see DATA_SOURCE)');


%%%%%%% range variables

id_range = netcdf.defVar(ncid,'range','nc_float',did_range);
netcdf.putAtt(ncid,id_range,'standard_name','RANGE');
netcdf.putAtt(ncid,id_range,'long_name','RANGE gate of the radar');
netcdf.putAtt(ncid,id_range,'units','m');
netcdf.putAtt(ncid,id_range,'valid_range', [min(data.range(:)), max(data.range(:))]);
netcdf.putAtt(ncid,id_range,'comment','Range from antenna to the center of each range gate');


%%%%%%%% chirp_seq_dependent variables

id_SeqAvg = netcdf.defVar(ncid,'seq_avg','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_SeqAvg,'standard_name','radar.operation.parameter.avg_chirps_per_chirp');
netcdf.putAtt(ncid,id_SeqAvg,'long_name','Number of averaged chirps in each chirp sequence');

id_SeqIntTime = netcdf.defVar(ncid,'seq_int_time','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_SeqIntTime,'standard_name','radar.operation.parameter.int_time_chirp');
netcdf.putAtt(ncid,id_SeqIntTime,'long_name','Integration time of each chirp sequence');
netcdf.putAtt(ncid,id_SeqIntTime,'units','seconds');

id_DoppLen = netcdf.defVar(ncid,'dopp_len','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_DoppLen,'standard_name','radar.operation.parameter.spec_samples_per_chirp');
netcdf.putAtt(ncid,id_DoppLen,'long_name','Number of samples in Dopppler spectra of each chirp sequence');
netcdf.putAtt(ncid,id_DoppLen,'comment','Needed to calculate the Doppler resolution: DoppRes = 2*nqv/dopp_len');

id_DoppMax = netcdf.defVar(ncid,'nqv','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_DoppMax,'standard_name','NYQUIST.VELOCITY');
netcdf.putAtt(ncid,id_DoppMax,'long_name',['Max. unambigious Doppler velocity for each chirp sequence; '...
                                          'Nyquist velocity per chirp sequence']);
netcdf.putAtt(ncid,id_DoppMax,'units','m s-1');
netcdf.putAtt(ncid,id_DoppMax,'comment','Needed to calculate the Doppler resolution: DoppRes = 2*nqv/dopp_len');

id_nAvg = netcdf.defVar(ncid,'n_avg','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_nAvg,'standard_name','radar.operation.parameter.no_spec_avg');
netcdf.putAtt(ncid,id_nAvg,'long_name','Number of spectra averaged for each chirp');
netcdf.putAtt(ncid,id_nAvg,'comment','n_avg = seq_avg/dopp_len')

id_range_offsets = netcdf.defVar(ncid,'range_offsets','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_range_offsets,'standard_name','radar.operation.parameter.range_index_chirp');
netcdf.putAtt(ncid,id_range_offsets,'long_name','Chirp sequence start index array in range array');
netcdf.putAtt(ncid,id_range_offsets,'comment','The command range(range_offsets) will give you the range where a new chirp sequence starts. range_offsets counts from 1 to n_levels.');


%%%%%%%% time dependend variables

id_time = netcdf.defVar(ncid,'time','nc_uint',did_time);
netcdf.putAtt(ncid,id_time,'standard_name','DATETIME');
netcdf.putAtt(ncid,id_time,'long_name','DATETIME in UTC');
netcdf.putAtt(ncid,id_time,'units','MJD2K');
netcdf.putAtt(ncid,id_time,'valid_range', [min(data.time(:)), max(data.time(:))]);
netcdf.putAtt(ncid,id_time,'comment','Time in sec since 2001.01.01. 00:00:00.');
if isfield(data, 'totsampchangelabel' )
    netcdf.putAtt(ncid,id_time, 'quality_flag', 'Dublicate time stamps found in lv0-file, the first occurrence of the dublicate time is removed')
end

id_sampleTms = netcdf.defVar(ncid,'sample_tms','nc_int',did_time);
netcdf.putAtt(ncid,id_sampleTms,'standarr_name','DATETIME.milliseconds');
netcdf.putAtt(ncid,id_sampleTms,'long_name','Milliseconds of sample');
netcdf.putAtt(ncid,id_sampleTms,'units','mu s');
netcdf.putAtt(ncid,id_sampleTms,'comment','To get the correct time the variable sample_tms must be added: time = time + sample_tms.');

id_RR = netcdf.defVar(ncid,'rr','nc_float',did_time);
netcdf.putAtt(ncid,id_RR,'standard_name','RAIN.RATE.SURFACE');
netcdf.putAtt(ncid,id_RR,'long_name','RAIN.RATE of meteo-station of the radar system');
netcdf.putAtt(ncid,id_RR,'units','mm h-1');
netcdf.putAtt(ncid,id_RR,'valid_range',[min(data.RR(:)), max(data.RR(:))]);
netcdf.putAtt(ncid,id_RR,'fill_value','NaNf');
netcdf.putAtt(ncid,id_RR,'rain.rate.surface_source','Vaisala weather station WXT520 or WXT530');

id_rh = netcdf.defVar(ncid,'rh','nc_float',did_time);
netcdf.putAtt(ncid,id_rh,'standard_name','HUMIDITY.RELATIVE.SURFACE');
netcdf.putAtt(ncid,id_rh,'long_name','Relative humidity of meteo-station');
netcdf.putAtt(ncid,id_rh,'units','%');
netcdf.putAtt(ncid,id_rh,'valid_range',[min(data.rh(:)), max(data.rh(:))]);
netcdf.putAtt(ncid,id_rh,'fill_value','NaNf');
netcdf.putAtt(ncid,id_rh,'humidity.relative.surface_source','Vaisala weather station WXT520 or WXT530');

id_T_env = netcdf.defVar(ncid,'ta','nc_float',did_time);
netcdf.putAtt(ncid,id_T_env,'standard_name','TEMPERATURE.SURFACE');
netcdf.putAtt(ncid,id_T_env,'long_name','TEMPERATURE.SURFACE of the environment measured by the meteo-station');
netcdf.putAtt(ncid,id_T_env,'units','K');
netcdf.putAtt(ncid,id_T_env,'valid_range',[min(data.T_env(:)), max(data.T_env(:))]);
netcdf.putAtt(ncid,id_T_env,'fill_value','NaNf');
netcdf.putAtt(ncid,id_T_env,'temperature.surface_source','Vaisala weather station WXT520 or WXT530');

id_pres = netcdf.defVar(ncid,'pa','nc_float',did_time);
netcdf.putAtt(ncid,id_pres,'standard_name','SURFACE.PRESSURE');
netcdf.putAtt(ncid,id_pres,'long_name','SURFACE.PRESSURE of the enviroment measured by the meteo-station');
netcdf.putAtt(ncid,id_pres,'units','hPa');
netcdf.putAtt(ncid,id_pres,'valid_range',[min(data.pres(:)), max(data.pres(:))]);
netcdf.putAtt(ncid,id_pres,'fill_value','NaNf');
netcdf.putAtt(ncid,id_pres,'surface.pressure_source','Vaisala weather station WXT520 or WXT530');

id_ff = netcdf.defVar(ncid,'wspeed','nc_float',did_time);
netcdf.putAtt(ncid,id_ff,'strandard_name','WIND.SPEED.SURFACE');
netcdf.putAtt(ncid,id_ff,'long_name','WIND.SPEED measured at about 1.5 m by the meteo-station');
netcdf.putAtt(ncid,id_ff,'units','m s-1');
netcdf.putAtt(ncid,id_ff,'valid_range',[min(data.ff(:)), max(data.ff(:))]);
netcdf.putAtt(ncid,id_ff,'fill_value','NaNf');
netcdf.putAtt(ncid,id_ff,'wind.speed.surface_source','Vaisala weather station WXT520 or WXT530');

id_fff = netcdf.defVar(ncid,'wdir','nc_float',did_time);
netcdf.putAtt(ncid,id_fff,'standard_name','WIND.DIRECTION.SURFACE');
netcdf.putAtt(ncid,id_fff,'long_name','WIND.DIRECTION measured at about 1.5 m by the meteo-station');
netcdf.putAtt(ncid,id_fff,'units','degrees');
netcdf.putAtt(ncid,id_fff,'valid_range',[min(data.fff(:)), max(data.fff(:))]);
netcdf.putAtt(ncid,id_fff,'fill_value','NaNf');
netcdf.putAtt(ncid,id_fff,'wind.direction.surface_source',['Vaisala weather '...
                           'station WXT520 or WXT530']);

id_Tb = netcdf.defVar(ncid,'tb','nc_float',did_time);
netcdf.putAtt(ncid,id_Tb,'standard_name','TEMPERATURE.BRIGHTNESS');
netcdf.putAtt(ncid,id_Tb,'long_name','Brightness temperature direct detection channel');
netcdf.putAtt(ncid,id_Tb,'units','K');
netcdf.putAtt(ncid,id_Tb,'valid_range',[min(data.Tb(:)), max(data.Tb(:))]);
netcdf.putAtt(ncid,id_Tb,'fill_value','NaNf');
netcdf.putAtt(ncid,id_Tb,'comment',['Brightness Temperature measurements from '...
                                    'the Passive 89-GHz Chanal of the radar']);

id_lwp = netcdf.defVar(ncid,'lwp','nc_float',did_time);
netcdf.putAtt(ncid,id_lwp,'standard_name','LIQUID.WATER.PATH');
netcdf.putAtt(ncid,id_lwp,'long_name','Liquid water path (lwp) calculated by RPG software');
netcdf.putAtt(ncid,id_lwp,'units','g m-2');
netcdf.putAtt(ncid,id_lwp,'valid_range',[min(data.lwp(:)), max(data.lwp(:))]);
netcdf.putAtt(ncid,id_lwp,'fill_value','NaNf');
netcdf.putAtt(ncid,id_lwp,'comment',['Liquid water path is calculated from '...
                                     'the tb measurement of the 89-GHz chanal. '...
                                     'The retrieval is developed by RPG and '...
                                     'based on a nural network approach']);

id_status = netcdf.defVar(ncid,'status','nc_float',did_time);
netcdf.putAtt(ncid,id_status,'standard_name','radar.operation.quality.flag.blower_status');
netcdf.putAtt(ncid,id_status,'long_name','blower status flag');
netcdf.putAtt(ncid,id_status,'comment',['The quality flag shows the following '...
                                        'status: 0/1 = heater on/off; '...
                                        '0/10 = blower on/off. The parameter '...
                                        'is recorded to check and monitor '...
                                        'the radar performance']);

id_TransPow = netcdf.defVar(ncid,'p_trans','nc_float',did_time);
netcdf.putAtt(ncid,id_TransPow,'standard_name','radar.operation.parameter.p_trans');
netcdf.putAtt(ncid,id_TransPow,'long_name','Transmitted power');
netcdf.putAtt(ncid,id_TransPow,'units','W');
netcdf.putAtt(ncid,id_TransPow,'valid_range',[min(data.TransPow(:)), max(data.TransPow(:))]);
netcdf.putAtt(ncid,id_TransPow,'fill_value','NaNf');
netcdf.putAtt(ncid,id_TransPow,'comment',['The transmitted power is is recorded '...
                                          'to check and monitor the radar performance '...
                                          'and measurement quality']);
                                 
id_T_trans = netcdf.defVar(ncid,'t_trans','nc_float',did_time);
netcdf.putAtt(ncid,id_T_trans,'standard_name','radar.operation.parameter.t_trans');
netcdf.putAtt(ncid,id_T_trans,'long_name','Transmitter temperature');
netcdf.putAtt(ncid,id_T_trans,'units','K');
netcdf.putAtt(ncid,id_T_trans,'valid_range',[min(data.T_trans(:)), max(data.T_trans(:))]);
netcdf.putAtt(ncid,id_T_trans,'fill_value','NaNf');
netcdf.putAtt(ncid,id_T_trans,'comment',['The transmitter temperature is recorded '...
                                         'to check and monitor the radar performance '...
                                          'and measurement quality']);
                                 
id_T_rec = netcdf.defVar(ncid,'t_rec','nc_float',did_time);
netcdf.putAtt(ncid,id_T_rec,'standard_name','radar.operation.parameter.t_rec');
netcdf.putAtt(ncid,id_T_rec,'long_name','Receiver temperature');
netcdf.putAtt(ncid,id_T_rec,'units','K');
netcdf.putAtt(ncid,id_T_rec,'valid_range',[min(data.T_rec(:)), max(data.T_rec(:))]);
netcdf.putAtt(ncid,id_T_rec,'fill_value','NaNf');
netcdf.putAtt(ncid,id_T_rec,'comment',['The receiver temperature is recorded '...
                                     'to check and monitor the radar performance '...
                                     'and measurement quality']);
                                 
id_T_pc = netcdf.defVar(ncid,'t_pc','nc_float',did_time);
netcdf.putAtt(ncid,id_T_pc,'standard_name','radar.operation.parameter.t_pc');
netcdf.putAtt(ncid,id_T_pc,'long_name','radar PC temperature');
netcdf.putAtt(ncid,id_T_pc,'units','K');
netcdf.putAtt(ncid,id_T_pc,'valid_range',[min(data.T_pc(:)), max(data.T_pc(:))]);
netcdf.putAtt(ncid,id_T_pc,'fill_value','NaNf');
netcdf.putAtt(ncid,id_T_pc,'comment',['The radar PC temperature is recorded '...
                                     'to check and monitor the radar performance '...
                                     'and measurement quality']);

id_QF = netcdf.defVar(ncid,'qf_rad','nc_byte',did_time);
netcdf.putAtt(ncid,id_QF,'standard_name','radar.operation.quality.flag.qf_radar');
netcdf.putAtt(ncid,id_QF,'long_name','Quality flag given by radar');
netcdf.putAtt(ncid,id_QF,'comment', ['To get the bit entries, one has to'...
                                     'convert the integer into a 4 bit binary. '...
                                     'bit4 = ADC saturation, bit3 = spectral '...
                                     'width too high, bit2 = no transmitter '...
                                     'power leveling. Note that in the above '...
                                     'convention holds: bit1 = 2^3, '...
                                     'bit2 = 2^2, bit3 = 2^1, bit4 = 2^0'])
                                 
%%%%%%%% multi-D variables

id_Ze = netcdf.defVar(ncid,'ze','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_Ze,'standard_name','RADAR.REFLECTIVITY.FACTOR_VV');
netcdf.putAtt(ncid,id_Ze,'long_name','Equivalent radar reflectivity factor Ze');
netcdf.putAtt(ncid,id_Ze,'units','mm6 m-3');
netcdf.putAtt(ncid,id_Ze,'valid_range',[min(data.Ze(:)), max(data.Ze(:))]);
netcdf.putAtt(ncid,id_Ze,'fill_value','NaNf');
netcdf.putAtt(ncid,id_Ze,'long_name',['The equivalent radar reflectivity '...
                          'factor Ze is obtained at vertical polarisation. '...
                          'If more polarisation states could be measured '...
                          'ze would be ze_vv']);
if isfield(data, 'Ze_label') % Ze corrected, adding note
    netcdf.putAtt(ncid,id_Ze,'comment',data.Ze_label);
    netcdf.putAtt(ncid,id_Ze,'corretion_dB',data.Ze_corr);
end

id_vm = netcdf.defVar(ncid,'vm','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_vm,'standard_name','DOPPLER.VELOCITY_MEAN');
netcdf.putAtt(ncid,id_vm,'long_name','Mean DOPPLER.VELOCITY');
netcdf.putAtt(ncid,id_vm,'units','m s-1');
netcdf.putAtt(ncid,id_vm,'valid_range',[min(data.vm(:)), max(data.vm(:))]);
netcdf.putAtt(ncid,id_vm,'fill_value','NaNf');
netcdf.putAtt(ncid,id_vm,'comment',['radial velocities of scatterers, negative '...
                                    'velocities indicate particles motion '...
                                    'towards the radar'])

id_sigma = netcdf.defVar(ncid,'sw','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_sigma,'standard_name','DOPPLER.SPECTRUM_WIDTH');
netcdf.putAtt(ncid,id_sigma,'long_name','Spectral width of Doppler velocity spectrum');
netcdf.putAtt(ncid,id_sigma,'units','m s-1');
netcdf.putAtt(ncid,id_sigma,'valid_range',[min(data.sigma(:)), max(data.sigma(:))]);
netcdf.putAtt(ncid,id_sigma,'fill_value','NaNf');
netcdf.putAtt(ncid,id_sigma,'comment',['Scectral width of Doppler spectrum at vertical '...
                                    'polarisation'])

id_skew = netcdf.defVar(ncid,'skew','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_skew,'standard_name','DOPPLER.SPECTRUM_SKEWNESS');
netcdf.putAtt(ncid,id_skew,'long_name','Doppler spectrum skewness');
netcdf.putAtt(ncid,id_skew,'valid_range',[min(data.skew(:)), max(data.skew(:))]);
netcdf.putAtt(ncid,id_skew,'fill_value','NaNf');
netcdf.putAtt(ncid,id_skew,'comment',['Doppler spectrum skewness of Doppler spectrum '...
                                    'at vertical polarisation'])
%Included by Bravo-Aranda, J.A. JABA
id_ldr = netcdf.defVar(ncid,'ldr','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_ldr,'standard_name','LINEAR.DEPOLARIZATION.RATIO');
netcdf.putAtt(ncid,id_ldr,'long_name','Linear depolarization ratio');
netcdf.putAtt(ncid,id_ldr,'unite','mm6 m-3');
netcdf.putAtt(ncid,id_ldr,'valid_range',[min(data.LDR(:)), max(data.LDR(:))]);
netcdf.putAtt(ncid,id_ldr,'fill values','NaNf');
netcdf.putAtt(ncid,id_ldr,'comment','L_dr = Ze_hv/Ze in [mm6 m-3]');
if data.DualPol == 2
    id_xcorr = netcdf.defVar(ncid,'rho_hv','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_xcorr,'standard_name','CRORRELATION.COEFFICIENT.');
    netcdf.putAtt(ncid,id_xcorr,'long_name','co-cross-channel correlation coefficient');
    netcdf.putAtt(ncid,id_xcorr,'valid_range',[min(data.xcorr(:)), max(data.xcorr(:))]);
    netcdf.putAtt(ncid,id_xcorr,'fill_value','NaNf');
    
    id_difphase = netcdf.defVar(ncid,'phi_dp','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_difphase,'standard_name','DIFFERENTIAL.PHASE');
    netcdf.putAtt(ncid,id_difphase,'long_name','co-cross-channel differential phase');   
    netcdf.putAtt(ncid,id_difphase,'unite','degree');   
    netcdf.putAtt(ncid,id_difphase,'valid_range',[min(data.difphase(:)), max(data.difphase(:))]); 
    netcdf.putAtt(ncid,id_difphase,'fill_value','NaNf'); 
end

id_QualFlag = netcdf.defVar(ncid,'qf_pro','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_QualFlag,'standard_name','processing.quality.flag.anti_alias');
netcdf.putAtt(ncid,id_QualFlag,'long_name','Quality flag, added in the additional data processing to alert for known issues');
netcdf.putAtt(ncid,id_QualFlag,'comment', ...
    ['This variable contains information on anything that might impact the quality ', ...
    'of the data at each pixel. Must be converted into three bit binary string. ', ...
	'If 0, i.e. dec2bin(QualityFlag,3) = 000, none of the included issues were ', ...
    'found. The definitions of each bit are given in the definition attribute.']);
netcdf.putAtt(ncid,id_QualFlag,'definition', ...
    ['If 2^0 bit is 1: this range gate is known to have aritifical spikes occurring ', ...
     'If 2^1 bit is 1: aircraft or other known flying non-meteorological object ', ...
     'If 2^2 bit is 1: wet-radome (was a problem for mirac-a for a time period '...
     'when coating missing)']);

id_Aliasmask = netcdf.defVar(ncid,'alias_mask','nc_byte',[did_range,did_time]);
netcdf.putAtt(ncid,id_Aliasmask,'standard_name','processing.quality.flag.alias_mask');
netcdf.putAtt(ncid,id_Aliasmask,'long_name','Mask array indicating in which bin dealiasing was applied');
netcdf.putAtt(ncid,id_Aliasmask,'comment',['The mask shows, if AnitAlias = 1, '...
                                'then dealiasing was applied by RPG software: '...
                                '0 = not applied; 1 = applied; '...
                                'If AntiAlias = 2, then dealiasing was applied '...
                                'in post-processing: 0 = no aliasing detected, '...
                                '1 = aliasing detected; if any bin equals 1 '...
                                '(while AntiAlias = 2) then the full column was dealiased.']);


%% ###################### initialize compression of all floats:
netcdf.defVarDeflate(ncid,id_RR,true,true,9);
netcdf.defVarDeflate(ncid,id_rh,true,true,9);
netcdf.defVarDeflate(ncid,id_T_env,true,true,9);
netcdf.defVarDeflate(ncid,id_pres,true,true,9);
netcdf.defVarDeflate(ncid,id_ff,true,true,9);
netcdf.defVarDeflate(ncid,id_fff,true,true,9);
netcdf.defVarDeflate(ncid,id_Tb,true,true,9);
netcdf.defVarDeflate(ncid,id_lwp,true,true,9);
netcdf.defVarDeflate(ncid,id_status,true,true,9);
netcdf.defVarDeflate(ncid,id_TransPow,true,true,9);
netcdf.defVarDeflate(ncid,id_T_trans,true,true,9);
netcdf.defVarDeflate(ncid,id_T_rec,true,true,9);
netcdf.defVarDeflate(ncid,id_T_pc,true,true,9);
netcdf.defVarDeflate(ncid,id_QF,true,true,9);
netcdf.defVarDeflate(ncid,id_Ze,true,true,9);
netcdf.defVarDeflate(ncid,id_vm,true,true,9);
netcdf.defVarDeflate(ncid,id_sigma,true,true,9);
netcdf.defVarDeflate(ncid,id_skew,true,true,9);
netcdf.defVarDeflate(ncid,id_QualFlag,true,true,9);
netcdf.defVarDeflate(ncid,id_Aliasmask,true,true,9);
if data.DualPol > 0
    netcdf.defVarDeflate(ncid,id_ldr,true,true,9); %JABA    
end 
    
if data.DualPol == 2
    netcdf.defVarDeflate(ncid,id_difphase,true,true,9); %JABA    
    netcdf.defVarDeflate(ncid,id_xcorr,true,true,9); %JABA  
end

netcdf.endDef(ncid);



%% ####################### put variables into file

% scalars
netcdf.putVar(ncid,id_AntiAlias,0,data.AntiAlias);
netcdf.putVar(ncid,id_freq,0,data.freq * 1e9);
netcdf.putVar(ncid,id_lon,0,data.Lon);
netcdf.putVar(ncid,id_lat,0,data.Lat);
netcdf.putVar(ncid,id_MSL,0,data.MSL);

% range dependet
netcdf.putVar(ncid,id_range,0,data.n_levels,data.range);


% chrip seq dependent variables
netcdf.putVar(ncid,id_range_offsets,0,data.no_chirp_seq,data.range_offsets);
netcdf.putVar(ncid,id_SeqAvg,0,data.no_chirp_seq,data.SeqAvg);
netcdf.putVar(ncid,id_SeqIntTime,0,data.no_chirp_seq,data.SeqIntTime);
netcdf.putVar(ncid,id_DoppMax,0,data.no_chirp_seq,data.DoppMax);
netcdf.putVar(ncid,id_DoppLen,0,data.no_chirp_seq,data.DoppLen);
netcdf.putVar(ncid,id_nAvg,0,data.no_chirp_seq,data.nAvg);


% time dependent variables
netcdf.putVar(ncid,id_time,0,data.totsamp,data.time);
netcdf.putVar(ncid,id_sampleTms,0,data.totsamp,data.sampleTms);
netcdf.putVar(ncid,id_RR,0,data.totsamp,data.RR);
netcdf.putVar(ncid,id_rh,0,data.totsamp,data.rh);
netcdf.putVar(ncid,id_T_env,0,data.totsamp,data.T_env);
netcdf.putVar(ncid,id_pres,0,data.totsamp,data.pres);
netcdf.putVar(ncid,id_ff,0,data.totsamp,data.ff * (1000/3600));
netcdf.putVar(ncid,id_fff,0,data.totsamp,data.fff);
netcdf.putVar(ncid,id_Tb,0,data.totsamp,data.Tb);
netcdf.putVar(ncid,id_lwp,0,data.totsamp,data.lwp);
netcdf.putVar(ncid,id_status,0,data.totsamp,data.status);
netcdf.putVar(ncid,id_TransPow,0,data.totsamp,data.TransPow);
netcdf.putVar(ncid,id_T_trans,0,data.totsamp,data.T_trans);
netcdf.putVar(ncid,id_T_rec,0,data.totsamp,data.T_rec);
netcdf.putVar(ncid,id_T_pc,0,data.totsamp,data.T_pc);
netcdf.putVar(ncid,id_QF,0,data.totsamp,data.QF);


% multidimensional variables
netcdf.putVar(ncid,id_Ze,[0,0],[data.n_levels,data.totsamp],data.Ze');
netcdf.putVar(ncid,id_vm,[0,0],[data.n_levels,data.totsamp],data.vm');
netcdf.putVar(ncid,id_sigma,[0,0],[data.n_levels,data.totsamp],data.sigma');
netcdf.putVar(ncid,id_skew,[0,0],[data.n_levels,data.totsamp],data.skew');
netcdf.putVar(ncid,id_QualFlag,[0,0],[data.n_levels,data.totsamp],data.QualFlag');
netcdf.putVar(ncid,id_Aliasmask,[0,0],[data.n_levels,data.totsamp],data.Aliasmask');

if data.DualPol > 0
    netcdf.putVar(ncid,id_ldr,[0,0],[data.n_levels,data.totsamp],data.LDR'); %JABA
end

if data.DualPol == 2
    disp('WARNING! skipping writing id_xcorr and difphase into files, variables not available (yet??)')
    netcdf.putVar(ncid,id_xcorr,[0,0],[data.n_levels,data.totsamp],data.xcorr'); %JABA
    netcdf.putVar(ncid,id_difphase,[0,0],[data.n_levels,data.totsamp],data.difphase'); %JABA
end

netcdf.close(ncid);

end
