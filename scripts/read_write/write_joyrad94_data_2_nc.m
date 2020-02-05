function write_joyrad94_data_2_nc(data, outfile, config)
% this function writes joyrad94 data into netcdf4

%################## Create a netCDF file.

ncid = netcdf.create(outfile,'NETCDF4'); 


%% ################# Define dimensions

did_time = netcdf.defDim(ncid,'time',data.totsamp);
did_range = netcdf.defDim(ncid,'range',data.n_levels);
did_vel = netcdf.defDim(ncid,'velocity',max(data.DoppLen));
did_no_seq = netcdf.defDim(ncid,'chirp_sequences',data.no_chirp_seq);
did_scalar = netcdf.defDim(ncid,'scalar',1);

if ne(data.T_altcount,0)
    did_T_range = netcdf.defDim(ncid,'T_range',data.T_altcount);
end

if ne(data.H_altcount,0)
    did_H_range = netcdf.defDim(ncid,'H_range',data.H_altcount);
end

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

%--- from here on the variables are not GEOMS formate conforme

id_filecode = netcdf.defVar(ncid,'filecode','nc_int',did_scalar);
netcdf.putAtt(ncid,id_filecode,'long_name','filecode/version number');

id_samples = netcdf.defVar(ncid,'samples','nc_int',did_scalar);
netcdf.putAtt(ncid,id_samples,'long_name','No. samples');

id_levels = netcdf.defVar(ncid,'n_levels','nc_int',did_scalar);
netcdf.putAtt(ncid,id_levels,'long_name','Number of range levels');

id_T_levels = netcdf.defVar(ncid,'n_T_levels','nc_int',did_scalar);
netcdf.putAtt(ncid,id_T_levels,'long_name','Number of range levels of temp. profile');

id_H_levels = netcdf.defVar(ncid,'n_H_levels','nc_int',did_scalar);
netcdf.putAtt(ncid,id_H_levels,'long_name','Number of range levels of hum. profile');

id_no_chirp_seq = netcdf.defVar(ncid,'no_chirp_seq','nc_int',did_scalar);
netcdf.putAtt(ncid,id_no_chirp_seq,'long_name','Number of chirp sequences');
netcdf.putAtt(ncid,id_no_chirp_seq,'comment',...
    'The radar can be programmed to run different resolution modes for different layers spanning several range gates');

id_CalInt = netcdf.defVar(ncid,'CalInt','nc_int',did_scalar);
netcdf.putAtt(ncid,id_CalInt,'long_name','Sample interval for automated zero calibrations');
netcdf.putAtt(ncid,id_CalInt,'units','no. samples');

id_AntSep = netcdf.defVar(ncid,'AntSep','nc_float',did_scalar);
netcdf.putAtt(ncid,id_AntSep,'long_name','Separation of antenna axis');
netcdf.putAtt(ncid,id_AntSep,'units','m');

id_AntDia = netcdf.defVar(ncid,'AntDia','nc_float',did_scalar);
netcdf.putAtt(ncid,id_AntDia,'long_name','Antenna diameter');
netcdf.putAtt(ncid,id_AntDia,'units','m');

id_AntG = netcdf.defVar(ncid,'AntG','nc_float',did_scalar);
netcdf.putAtt(ncid,id_AntG,'long_name','Antenna gain');
netcdf.putAtt(ncid,id_AntG,'units','[]');

id_HPBW = netcdf.defVar(ncid,'HPBW','nc_float',did_scalar);
netcdf.putAtt(ncid,id_HPBW,'long_name','Cassegrain antenna HPBW');
netcdf.putAtt(ncid,id_HPBW,'units','degrees');

id_SampDur = netcdf.defVar(ncid,'SampDur','nc_float',did_scalar);
netcdf.putAtt(ncid,id_SampDur,'long_name','Full sample duration');
netcdf.putAtt(ncid,id_SampDur,'units','seconds');
netcdf.putAtt(ncid,id_SampDur,'comment','This is not the integration time (see SeqIntTime), but the sum of the chirp sequences integration times.')

id_C = netcdf.defVar(ncid,'C','nc_float',did_scalar);
netcdf.putAtt(ncid,id_C,'long_name','Radar constant (see manual eq. 2.1.5)');

id_pol = netcdf.defVar(ncid,'DualPol','nc_byte',did_scalar);
netcdf.putAtt(ncid,id_pol,'long_name','Polarisation. DualPol = 1: Polarimetric radar. DualPol = 0: Only vertical polarization available.');

if isfield(data, 'CompEna')
    id_compress = netcdf.defVar(ncid,'CompEna','nc_byte',did_scalar);
    netcdf.putAtt(ncid,id_compress,'long_name','If CompEna = 1, then compression was enabled in RPG software. Noise level from spectra was removed.');
end

id_cal_mom = netcdf.defVar(ncid,'cal_mom','nc_byte',did_scalar);
netcdf.putAtt(ncid,id_cal_mom,'long_name','Integer indicating how moments were calculated.');
netcdf.putAtt(ncid,id_cal_mom,'comment',...
    'Set in config file. 2 = moments were calculated from spectra. 3 = moments are taken from RPG processing sofware (lv1).');


%%%%%%% range variables

id_range = netcdf.defVar(ncid,'range','nc_float',did_range);
netcdf.putAtt(ncid,id_range,'standard_name','RANGE');
netcdf.putAtt(ncid,id_range,'long_name','RANGE gate of the radar');
netcdf.putAtt(ncid,id_range,'units','m');
netcdf.putAtt(ncid,id_range,'valid_range', [min(data.range(:)), max(data.range(:))]);
netcdf.putAtt(ncid,id_range,'comment','Range from antenna to the center of each range gate');

%--- from here on the variables are not GEOMS formate conforme

id_Fr = netcdf.defVar(ncid,'Fr','nc_int',did_range);
netcdf.putAtt(ncid,id_Fr,'long_name','Range factors, see manual eq. (2.5.6)');

if ne(data.T_altcount,0)
    id_T_levels = netcdf.defVar(ncid,'T_levels','nc_float',did_T_range);
    netcdf.putAtt(ncid,id_T_levels,'long_name','No. levels of temperature profile');
end

if ne(data.H_altcount,0)
    id_H_levels = netcdf.defVar(ncid,'H_levels','nc_float',did_H_range);
    netcdf.putAtt(ncid,id_H_levels,'long_name','No. levels of humidity profiles');
end


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

id_dr = netcdf.defVar(ncid,'dz_chirp','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_dr,'lstandard_name','radar.operation.parameter.dz_chirp');
netcdf.putAtt(ncid,id_dr,'long_name','Range resolution for chirp sequences');
netcdf.putAtt(ncid,id_dr,'units','m');


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

id_vol = netcdf.defVar(ncid,'vol','nc_float',did_time);
netcdf.putAtt(ncid,id_vol,'standard_name','radar.operation.parameter.vol');
netcdf.putAtt(ncid,id_vol,'long_name','Voltage direct detection channel');
netcdf.putAtt(ncid,id_vol,'units','V');

id_powIF = netcdf.defVar(ncid,'pow_if','nc_float',did_time);
netcdf.putAtt(ncid,id_powIF,'standard_name','radar.operation.parameter.pow_if');
netcdf.putAtt(ncid,id_powIF,'long_name','IF power at ADC');
netcdf.putAtt(ncid,id_powIF,'units','microWatts');

id_ele = netcdf.defVar(ncid,'ele','nc_float',did_time);
netcdf.putAtt(ncid,id_ele,'standard_name','radar.operation.parameter.system_elevation');
netcdf.putAtt(ncid,id_ele,'long_name','Sensor elevation angle');
netcdf.putAtt(ncid,id_ele,'units','degrees');

id_az = netcdf.defVar(ncid,'azi','nc_float',did_time);
netcdf.putAtt(ncid,id_az,'standard_name','radar.operation.parameter.system_azimuth');
netcdf.putAtt(ncid,id_az,'long_name','Sensor azimuth angle');
netcdf.putAtt(ncid,id_az,'units','degrees');

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

id_incel = netcdf.defVar(ncid,'incl_el','nc_float',did_time);
netcdf.putAtt(ncid,id_incel,'long_name','Inclination angle; elevation');
netcdf.putAtt(ncid,id_incel,'units','min');

id_incea = netcdf.defVar(ncid,'incl_ea','nc_float',did_time);
netcdf.putAtt(ncid,id_incea,'long_name','Inclination angle; angle of elevation axis');
netcdf.putAtt(ncid,id_incea,'units','min');


%%%%%%%% multi-D variables

id_vel = netcdf.defVar(ncid,'chirp_vel','nc_float',[did_vel,did_no_seq]);
netcdf.putAtt(ncid,id_vel,'long_name','Velocity arrays for each chirp sequence');
netcdf.putAtt(ncid,id_vel,'units','m/s');
netcdf.putAtt(ncid,id_vel,'comment',...
    ['If spectra has been dealiased, use MinVel to get the correct velocity array for each spectra.' ...
    'The velocity array is asymmetric, i.e. the absolute values of maximum and minumum velocities are not equal. ',...
    'Since the spectrum at -v_nyquist and +v_nyquist is the same, the entry at +v_nyquist was cut.']);


id_mask = netcdf.defVar(ncid,'mask','nc_byte',[did_range,did_time]);
netcdf.putAtt(ncid,id_mask,'long_name','Mask array of occupied range cells: 0 = not occupied; 1 = occupied');
netcdf.putAtt(ncid,id_mask,'comment',['Mask(time,range) is set to 0 when signal does not exceed threshold set in the software. ',...
    'The threshold is a multiple of the noise standard deviation that is set in the measurement definition file and defines ',...
    'the linear sensitivity limit (see variable SLv). Doppler spectrum is only stored when set to 1. ',...
    'Note that data produced with software version 1 (i.e. before 7 Feb 2017) used a different threshold, ',...
    'i.e. the standard deviation for each chrip sequence. ',...
    'Moreover, the the latter was calculated wrongly by the software version 1, which lead to the rejection of actual valid samples.']);

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
                            
id_AliasStatus = netcdf.defVar(ncid,'AliasStatus','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_AliasStatus,'long_name','Flags indicating quality of dealiasing.');
netcdf.putAtt(ncid,id_AliasStatus,'comment',['Must be converted into four bit binary string. If 0, i.e. dec2bin(AliasStatus) = 0000 then dealiasing was successful. ',...
    'If 2^0 bit is 1, there was no initial guess velocity for this bin. ',...
    'If 2^1 bit is 1, the spectrum where the main peak was assumed to occur could not be concatenated properly since chirp sequence boundaries were reached. ',...
    'If 2^2 bit is 1, there is still significant signal close to the nyquist limit. dealiasing might not have been performed properly. ',...
    'If 2^3 bit is 1, the mean of the mean Doppler velocity in differs by more than 5 m/s from the neighbouring bins mean value. ',...
    'The larger the decimal number, the less reliable the data. ',...
    'This flag is only available if dealiasing was applied in process_joyrad94_data.m. For RPG dealising there is no quality indicator.']);

id_MinVel = netcdf.defVar(ncid,'MinVel','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_MinVel,'long_name','Minimum velocity in velocity array');
netcdf.putAtt(ncid,id_MinVel,'comment',['If dealising was applied the array indicates the velocity of the first bin in spec. ',...
    'The correct velocity array can be obtained by: vel_true = velocity + Minvel(i,j) - velocity(1)']);
netcdf.putAtt(ncid,id_MinVel,'units','m/s');

id_MinVel_Correction = netcdf.defVar(ncid,'MinVel_Correction','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_MinVel_Correction,'long_name','Minimum velocity correction');
netcdf.putAtt(ncid,id_MinVel_Correction,'comment',['If dealising was applied the array indicates by how much MinVel '...
    'was corrected by the final quality check (output from dealias_spectra_vm_column_quality_check.m). '...
    'Subtracting this value from MinVel provides the offset before the final quality check (output from dealias_spectra.m).']);
netcdf.putAtt(ncid,id_MinVel_Correction,'units','m/s');

%--- radar moments ---

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

id_kurt = netcdf.defVar(ncid,'kurt','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_kurt,'standard_name','RADAR.DOPPLER.SPECTRUM_KURTOSIS');
netcdf.putAtt(ncid,id_kurt,'long_name','Doppler spectrum Kurtosis');
netcdf.putAtt(ncid,id_kurt,'valid_range',[min(data.kurt(:)), max(data.kurt(:))]);
netcdf.putAtt(ncid,id_kurt,'fill_value','NaNf');
netcdf.putAtt(ncid,id_kurt,'comment',['Doppler spectrum Kurtosis of Doppler spectrum '...
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
    %netcdf.putAtt(ncid,id_xcorr,'valid_range',[min(data.xcorr(:)), max(data.xcorr(:))]);
    netcdf.putAtt(ncid,id_xcorr,'fill_value','NaNf');
    
    id_difphase = netcdf.defVar(ncid,'phi_dp','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_difphase,'standard_name','DIFFERENTIAL.PHASE');
    netcdf.putAtt(ncid,id_difphase,'long_name','co-cross-channel differential phase');   
    netcdf.putAtt(ncid,id_difphase,'unite','degree');   
    %netcdf.putAtt(ncid,id_difphase,'valid_range',[min(data.difphase(:)), max(data.difphase(:))]); 
    netcdf.putAtt(ncid,id_difphase,'fill_value','NaNf');
end

id_PNv = netcdf.defVar(ncid,'PNv','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_PNv,'long_name','Total IF power in vertical polarization measured at ADC input');
netcdf.putAtt(ncid,id_PNv,'units','W');

id_SLv = netcdf.defVar(ncid,'SLv','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_SLv,'long_name','Linear sensitivity limit for vertical polarisation');
netcdf.putAtt(ncid,id_SLv,'units','mm6 mm-3');

id_spec = netcdf.defVar(ncid,'sze','nc_float',[did_vel,did_range,did_time]);
netcdf.putAtt(ncid,id_spec,'standard_name','RADAR.DOPPLER.SPECTRUM_VV');
netcdf.putAtt(ncid,id_spec,'long_name','Doppler spectrum vertical polarization');
netcdf.putAtt(ncid,id_spec,'units','mm6 mm-3');
netcdf.putAtt(ncid,id_spec,'valid_range',[min(data.spec), max(data.spec)]);
netcdf.putAtt(ncid,id_spec,'fill_value','NaNf');
netcdf.putAtt(ncid,id_spec,'comment',['This is the normalized Doppler spectrum. '...
                           'The integral of the spectra minus Ze gives the noise level, ',...
                           'i.e. integral(spec)-Ze = noise. The velocity array is asymmetric, '...
                           'i.e. the absolute values of maximum and minumum velocities are not equal. ',...
                           'Since the spectrum at -v_nyquist and +v_nyquist is the same, '...
                           'the entry at +v_nyquist was cut.']);

id_VNoisePow_mean = netcdf.defVar(ncid,'v_noise_pow_mean','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_VNoisePow_mean,'standard_name','radar.doppler.spectrum_noise.power.mean_vv');
netcdf.putAtt(ncid,id_VNoisePow_mean,'long_name','Bin Doppler spectrum mean noise power in vertical polarization');
netcdf.putAtt(ncid,id_VNoisePow_mean,'units','mm6 m-3');
netcdf.putAtt(ncid,id_VNoisePow_mean,'valid_range',[min(data.VNoisePow_mean), max(data.VNoisePow_mean)]);
netcdf.putAtt(ncid,id_VNoisePow_mean,'fill_value','NaNf');
netcdf.putAtt(ncid,id_VNoisePow_mean,'comment','If cal_mom == 3, then bin noise power is only calculated if CompEna == 1.');

id_VNoisePow_peak = netcdf.defVar(ncid,'v_noise_pow_peak','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_VNoisePow_peak,'standard_name','radar.doppler.spectrum_noise.power.peak_vv');
netcdf.putAtt(ncid,id_VNoisePow_peak,'long_name','Bin Doppler spectrum peak noise power in vertical polarization');
netcdf.putAtt(ncid,id_VNoisePow_peak,'units','mm6 m-3');
netcdf.putAtt(ncid,id_VNoisePow_peak,'valid_range',[min(data.VNoisePow_peak), max(data.VNoisePow_peak)]);
netcdf.putAtt(ncid,id_VNoisePow_peak,'fill_value','NaNf');
netcdf.putAtt(ncid,id_VNoisePow_peak,'comment','If cal_mom == 3, then bin noise power is only calculated if CompEna == 1.');

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




if data.DualPol > 0
    
    id_PNh = netcdf.defVar(ncid,'if_pow_hh','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_PNh,'long_name','Total IF power in horizontal polarization measured at ADC input');
    netcdf.putAtt(ncid,id_PNh,'units','W');
    
    id_SLh = netcdf.defVar(ncid,'sen_lim_vv','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_SLh,'long_name','Linear sensitivity limit for vertical polarisation');
    netcdf.putAtt(ncid,id_SLh,'units','mm6 mm-3');
    
    id_spec_h = netcdf.defVar(ncid,'sze_hh','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_spec_h,'long_name','Doppler spectrum horizontal polarization');
    netcdf.putAtt(ncid,id_spec_h,'units','mm6 mm-3');
    
    id_spec_covRe = netcdf.defVar(ncid,'scov_re','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_spec_covRe,'long_name','Real part of covariance spectrum');
    netcdf.putAtt(ncid,id_spec_covRe,'units','mm6 mm-3');
    
    id_spec_covIm = netcdf.defVar(ncid,'scov_im','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_spec_covIm,'long_name','Real part of covariance spectrum');
    netcdf.putAtt(ncid,id_spec_covIm,'units','mm6 mm-3');

end % if data.DualPol > 0


if data.CompEna == 2 && data.DualPol == 2
    
    id_d_spec = netcdf.defVar(ncid,'sz_dr','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_d_spec,'long_name','Differntial spectral reflectivity');
    netcdf.putAtt(ncid,id_d_spec,'units','mm6 mm-3');
    netcdf.putAtt(ncid,id_d_spec,'valid_range',[min(data.d_spec), max(data.d_spec)]);
    netcdf.putAtt(ncid,id_d_spec,'fill_value','NaNf');
    
    id_CorrCoeff = netcdf.defVar(ncid,'srho_hv','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_CorrCoeff,'long_name','Spectral correlation coefficient');
    netcdf.putAtt(ncid,id_CorrCoeff,'valid_range',[min(data.CorrCoeff), max(data.CorrCoeff)]);
    netcdf.putAtt(ncid,id_CorrCoeff,'fill_value','NaNf');
    
    id_DiffPh = netcdf.defVar(ncid,'sphi_dp','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_DiffPh,'long_name','Spectral differential phase');
    netcdf.putAtt(ncid,id_DiffPh,'units','degree');
    netcdf.putAtt(ncid,id_DiffPh,'valid_range',[min(data.DiffPh), max(data.DiffPh)]);
    netcdf.putAtt(ncid,id_DiffPh,'fill_value','NaNf');
    
end % if data.CompEna == 2 && data.DualPol > 0


if data.DualPol == 2 && data.CompEna == 2
    
    id_SLDR = netcdf.defVar(ncid,'slt_ldr','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_SLDR,'long_name','Spectral slanted LDR');
    netcdf.putAtt(ncid,id_SLDR,'units','mm6 mm-3');
    netcdf.putAtt(ncid,id_SLDR,'valid_range',[min(data.SLDR), max(data.SLDR)]);
    netcdf.putAtt(ncid,id_SLDR,'fill_value','NaNf');
    
    id_SCorrCoeff = netcdf.defVar(ncid,'sslt_ldr','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_SCorrCoeff,'long_name','Spectral slanted corellation coefficient');
    netcdf.putAtt(ncid,id_SCorrCoeff,'valid_range',[min(data.SCorrCoeff), max(data.SCorrCoeff)]);
    netcdf.putAtt(ncid,id_SCorrCoeff,'fill_value','NaNf');
    
    if data.CompEna == 2
        
        id_KDP = netcdf.defVar(ncid,'kdp','nc_float',[did_range,did_time]);
        netcdf.putAtt(ncid,id_KDP,'long_name','Specific differential phase shift');
        netcdf.putAtt(ncid,id_KDP,'units','degree km-1');
        netcdf.putAtt(ncid,id_KDP,'valid_range',[min(data.KDP), max(data.KDP)]);
        netcdf.putAtt(ncid,id_KDP,'fill_value','NaNf');
        
        id_DiffAtt = netcdf.defVar(ncid,'att_dif','nc_float',[did_range,did_time]);
        netcdf.putAtt(ncid,id_DiffAtt,'long_name','Differential attenuation');
        netcdf.putAtt(ncid,id_DiffAtt,'units','mm6 mm-3 km-1');    
        netcdf.putAtt(ncid,id_DiffAtt,'valid_range',[min(data.DiffAtt), max(data.DiffAtt)]);
        netcdf.putAtt(ncid,id_DiffAtt,'fill_value','NaNf');        
        
    end % if data.CompEna == 2
    
end % if data.DualPol == 2 && data.CompEna == 2

if data.CompEna > 0 && data.DualPol > 0
         id_HNoisePow_mean = netcdf.defVar(ncid,'h_noise_pow_mean','nc_float',[did_range,did_time]);
         netcdf.putAtt(ncid,id_HNoisePow_mean,'standard_name','radar.doppler.spectrum_noise.power.mean_hh');         
         netcdf.putAtt(ncid,id_HNoisePow_mean,'long_name','Bin Doppler spectrum mean noise power in horizontal polarization');
         netcdf.putAtt(ncid,id_HNoisePow_mean,'units','mm6 mm-3');
         netcdf.putAtt(ncid,id_HNoisePow_mean,'comment','If cal_mom == 3, then bin noise power is only calculated if CompEna == 1.');
         
         id_HNoisePow_peak = netcdf.defVar(ncid,'h_noise_pow_peak','nc_float',[did_range,did_time]);
         netcdf.putAtt(ncid,id_HNoisePow_peak,'standard_name','radar.doppler.spectrum_noise.power.peak_hh');                  
         netcdf.putAtt(ncid,id_HNoisePow_peak,'long_name','Bin Doppler spectrum peak noise power in horizontal polarization');
         netcdf.putAtt(ncid,id_HNoisePow_peak,'units','mm6 mm-3');
         netcdf.putAtt(ncid,id_HNoisePow_peak,'comment','If cal_mom == 3, then bin noise power is only calculated if CompEna == 1.');
end % data.CompEna > 0


%% ###################### initialize copression
netcdf.defVarDeflate(ncid,id_RR,true,true,9);
netcdf.defVarDeflate(ncid,id_rh,true,true,9);
netcdf.defVarDeflate(ncid,id_T_env,true,true,9);
netcdf.defVarDeflate(ncid,id_pres,true,true,9);
netcdf.defVarDeflate(ncid,id_ff,true,true,9);
netcdf.defVarDeflate(ncid,id_fff,true,true,9);
netcdf.defVarDeflate(ncid,id_vol,true,true,9);
netcdf.defVarDeflate(ncid,id_Tb,true,true,9);
netcdf.defVarDeflate(ncid,id_lwp,true,true,9);
netcdf.defVarDeflate(ncid,id_powIF,true,true,9);
netcdf.defVarDeflate(ncid,id_ele,true,true,9);
netcdf.defVarDeflate(ncid,id_az,true,true,9);
netcdf.defVarDeflate(ncid,id_status,true,true,9);
netcdf.defVarDeflate(ncid,id_TransPow,true,true,9);
netcdf.defVarDeflate(ncid,id_T_trans,true,true,9);
netcdf.defVarDeflate(ncid,id_T_rec,true,true,9);
netcdf.defVarDeflate(ncid,id_T_pc,true,true,9);
netcdf.defVarDeflate(ncid,id_QF,true,true,9);
netcdf.defVarDeflate(ncid,id_incel,true,true,9);
netcdf.defVarDeflate(ncid,id_incea,true,true,9);
netcdf.defVarDeflate(ncid,id_mask,true,true,9);
netcdf.defVarDeflate(ncid,id_Aliasmask,true,true,9);
netcdf.defVarDeflate(ncid,id_AliasStatus,true,true,9);
netcdf.defVarDeflate(ncid,id_MinVel,true,true,9);
netcdf.defVarDeflate(ncid,id_Ze,true,true,9);
netcdf.defVarDeflate(ncid,id_vm,true,true,9);
netcdf.defVarDeflate(ncid,id_sigma,true,true,9);
netcdf.defVarDeflate(ncid,id_skew,true,true,9);
netcdf.defVarDeflate(ncid,id_kurt,true,true,9);
netcdf.defVarDeflate(ncid,id_PNv,true,true,9);
netcdf.defVarDeflate(ncid,id_SLv,true,true,9);
netcdf.defVarDeflate(ncid,id_spec,true,true,9);
netcdf.defVarDeflate(ncid,id_VNoisePow_mean,true,true,9);
netcdf.defVarDeflate(ncid,id_VNoisePow_peak,true,true,9);
netcdf.defVarDeflate(ncid,id_QualFlag,true,true,9);

if data.DualPol > 0
    netcdf.defVarDeflate(ncid,id_PNh,true,true,9);
    netcdf.defVarDeflate(ncid,id_SLh,true,true,9);
    netcdf.defVarDeflate(ncid,id_spec_h,true,true,9);
    netcdf.defVarDeflate(ncid,id_spec_covRe,true,true,9);
    netcdf.defVarDeflate(ncid,id_spec_covIm,true,true,9);
    netcdf.defVarDeflate(ncid,id_ldr,true,true,9); %JABA    
    if data.DualPol == 2
        netcdf.defVarDeflate(ncid,id_difphase,true,true,9); %JABA    
        netcdf.defVarDeflate(ncid,id_xcorr,true,true,9); %JABA    
    end
end

if data.CompEna == 2 && data.DualPol == 2
    netcdf.defVarDeflate(ncid,id_d_spec,true,true,9);
    netcdf.defVarDeflate(ncid,id_CorrCoeff,true,true,9);
    netcdf.defVarDeflate(ncid,id_DiffPh,true,true,9);
end

if data.DualPol == 2 && data.CompEna == 2
    netcdf.defVarDeflate(ncid,id_SLDR,true,true,9);
    netcdf.defVarDeflate(ncid,id_SCorrCoeff,true,true,9);
    if data.CompEna == 2
        netcdf.defVarDeflate(ncid,id_KDP,true,true,9);
        netcdf.defVarDeflate(ncid,id_DiffAtt,true,true,9);
    end
    
end

if data.CompEna > 0 && data.DualPol > 0
    netcdf.defVarDeflate(ncid,id_HNoisePow_mean,true,true,9);
    netcdf.defVarDeflate(ncid,id_HNoisePow_peak,true,true,9);
end    

netcdf.endDef(ncid);





%% ####################### put variables into file

% scalars
netcdf.putVar(ncid,id_filecode,0,data.filecode);
netcdf.putVar(ncid,id_samples,0,data.totsamp);
netcdf.putVar(ncid,id_levels,0,data.n_levels);
netcdf.putVar(ncid,id_no_chirp_seq,0,data.no_chirp_seq);
netcdf.putVar(ncid,id_T_levels,0,data.T_altcount);
netcdf.putVar(ncid,id_H_levels,0,data.H_altcount);
netcdf.putVar(ncid,id_CalInt,0,data.CalInt);
netcdf.putVar(ncid,id_AntSep,0,data.AntSep);
netcdf.putVar(ncid,id_AntDia,0,data.AntDia);
netcdf.putVar(ncid,id_AntG,0,data.AntG);
netcdf.putVar(ncid,id_HPBW,0,data.HPBW);
netcdf.putVar(ncid,id_SampDur,0,data.SampDur);
netcdf.putVar(ncid,id_C,0,data.C);
netcdf.putVar(ncid,id_pol,0,data.DualPol);
if isfield(data, 'CompEna')
    netcdf.putVar(ncid,id_compress,0,data.CompEna);
end
netcdf.putVar(ncid,id_AntiAlias,0,data.AntiAlias);
netcdf.putVar(ncid,id_cal_mom,0,data.cal_mom);
netcdf.putVar(ncid,id_freq,0,data.freq);
netcdf.putVar(ncid,id_lon,0,data.Lon);
netcdf.putVar(ncid,id_lat,0,data.Lat);
netcdf.putVar(ncid,id_MSL,0,data.MSL);

% range dependet
netcdf.putVar(ncid,id_range,0,data.n_levels,data.range);
netcdf.putVar(ncid,id_Fr,0,data.n_levels,data.Fr);
if ne(data.T_altcount,0)
    netcdf.putVar(ncid,id_T_levels,0,data.T_altcount,data.T_alt);
end
if ne(data.H_altcount,0)
    netcdf.putVar(ncid,id_H_levels,0,data.H_altcount,data.H_alt);
end


% chrip seq dependent variables
netcdf.putVar(ncid,id_range_offsets,0,data.no_chirp_seq,data.range_offsets);
netcdf.putVar(ncid,id_dr,0,data.no_chirp_seq,data.dr);
netcdf.putVar(ncid,id_DoppMax,0,data.no_chirp_seq,data.DoppMax);
netcdf.putVar(ncid,id_DoppLen,0,data.no_chirp_seq,data.DoppLen);
netcdf.putVar(ncid,id_SeqAvg,0,data.no_chirp_seq,data.SeqAvg);
netcdf.putVar(ncid,id_SeqIntTime,0,data.no_chirp_seq,data.SeqIntTime);
netcdf.putVar(ncid,id_nAvg,0,data.no_chirp_seq,data.nAvg);


% time dependent variables
netcdf.putVar(ncid,id_time,0,data.totsamp,data.time);
netcdf.putVar(ncid,id_sampleTms,0,data.totsamp,data.sampleTms);
netcdf.putVar(ncid,id_RR,0,data.totsamp,data.RR);
netcdf.putVar(ncid,id_rh,0,data.totsamp,data.rh);
netcdf.putVar(ncid,id_T_env,0,data.totsamp,data.T_env);
netcdf.putVar(ncid,id_pres,0,data.totsamp,data.pres);
netcdf.putVar(ncid,id_ff,0,data.totsamp,data.ff);
netcdf.putVar(ncid,id_fff,0,data.totsamp,data.fff);
netcdf.putVar(ncid,id_vol,0,data.totsamp,data.u);
netcdf.putVar(ncid,id_Tb,0,data.totsamp,data.Tb);
netcdf.putVar(ncid,id_lwp,0,data.totsamp,data.lwp);
netcdf.putVar(ncid,id_powIF,0,data.totsamp,data.powIF);
netcdf.putVar(ncid,id_ele,0,data.totsamp,data.ele);
netcdf.putVar(ncid,id_az,0,data.totsamp,data.az);
netcdf.putVar(ncid,id_status,0,data.totsamp,data.status);
netcdf.putVar(ncid,id_TransPow,0,data.totsamp,data.TransPow);
netcdf.putVar(ncid,id_T_trans,0,data.totsamp,data.T_trans);
netcdf.putVar(ncid,id_T_rec,0,data.totsamp,data.T_rec);
netcdf.putVar(ncid,id_T_pc,0,data.totsamp,data.T_pc);
netcdf.putVar(ncid,id_QF,0,data.totsamp,data.QF);
netcdf.putVar(ncid,id_incel,0,data.totsamp,data.reserved(:,2));
netcdf.putVar(ncid,id_incea,0,data.totsamp,data.reserved(:,3));

% multidimensional variables
netcdf.putVar(ncid,id_vel,[0,0],[max(data.DoppLen),data.no_chirp_seq],data.velocity');
netcdf.putVar(ncid,id_mask,[0,0],[data.n_levels,data.totsamp],data.mask');
netcdf.putVar(ncid,id_Aliasmask,[0,0],[data.n_levels,data.totsamp],data.Aliasmask');
netcdf.putVar(ncid,id_AliasStatus,[0,0],[data.n_levels,data.totsamp],data.AliasStatus');
netcdf.putVar(ncid,id_MinVel,[0,0],[data.n_levels,data.totsamp],data.MinVel');
netcdf.putVar(ncid,id_MinVel_Correction,[0,0],[data.n_levels,data.totsamp],data.MinVel_Correction');
netcdf.putVar(ncid,id_Ze,[0,0],[data.n_levels,data.totsamp],data.Ze');
netcdf.putVar(ncid,id_vm,[0,0],[data.n_levels,data.totsamp],data.vm');
netcdf.putVar(ncid,id_sigma,[0,0],[data.n_levels,data.totsamp],data.sigma');
netcdf.putVar(ncid,id_skew,[0,0],[data.n_levels,data.totsamp],data.skew');
netcdf.putVar(ncid,id_kurt,[0,0],[data.n_levels,data.totsamp],data.kurt');
netcdf.putVar(ncid,id_PNv,[0,0],[data.n_levels,data.totsamp],data.PNv');
netcdf.putVar(ncid,id_SLv,[0,0],[data.n_levels,data.totsamp],data.SLv');
netcdf.putVar(ncid,id_spec,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.spec,[3,2,1]));
netcdf.putVar(ncid,id_VNoisePow_mean,[0,0],[data.n_levels,data.totsamp],data.VNoisePow_mean');
netcdf.putVar(ncid,id_VNoisePow_peak,[0,0],[data.n_levels,data.totsamp],data.VNoisePow_peak');
netcdf.putVar(ncid,id_QualFlag,[0,0],[data.n_levels,data.totsamp],data.QualFlag');



if data.DualPol > 0
    netcdf.putVar(ncid,id_PNh,[0,0],[data.n_levels,data.totsamp],data.PNh');
    netcdf.putVar(ncid,id_SLh,[0,0],[data.n_levels,data.totsamp],data.SLh');
    netcdf.putVar(ncid,id_spec_h,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.spec_hv,[3,2,1]));
    netcdf.putVar(ncid,id_spec_covRe,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.spec_covRe,[3,2,1]));
    netcdf.putVar(ncid,id_spec_covIm,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.spec_covIm,[3,2,1]));
    netcdf.putVar(ncid,id_ldr,[0,0],[data.n_levels,data.totsamp],data.LDR'); %JABA
    if data.DualPol == 2
        disp('WARNING! skipping writing id_xcorr and difphase into files, variables not available (yet??)')
%         netcdf.putVar(ncid,id_xcorr,[0,0],[data.n_levels,data.totsamp],data.xcorr'); %JABA
%         netcdf.putVar(ncid,id_difphase,[0,0],[data.n_levels,data.totsamp],data.difphase'); %JABA
    end
end

if data.CompEna == 2 && data.DualPol == 2
    netcdf.putVar(ncid,id_d_spec,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.d_spec,[3,2,1]));
    netcdf.putVar(ncid,id_CorrCoeff,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.CorrCoeff,[3,2,1]));
    netcdf.putVar(ncid,id_DiffPh,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.DiffPh,[3,2,1]));
end

if data.DualPol == 2 && data.CompEna == 2
    netcdf.putVar(ncid,id_SLDR,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.SLDR,[3,2,1]));
    netcdf.putVar(ncid,id_SCorrCoeff,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.SCorrCoeff,[3,2,1]));
    if data.CompEna == 2
        netcdf.putVar(ncid,id_KDP,[0,0],[data.n_levels,data.totsamp],data.KDP');
        netcdf.putVar(ncid,id_DiffAtt,[0,0],[data.n_levels,data.totsamp],data.DiffAtt');
    end
end

if data.CompEna > 0 && data.DualPol > 0
    netcdf.putVar(ncid,id_HNoisePow_mean,[0,0],[data.n_levels,data.totsamp],data.HNoisePow_mean');
    netcdf.putVar(ncid,id_HNoisePow_peak,[0,0],[data.n_levels,data.totsamp],data.HNoisePow_peak');
end

netcdf.close(ncid);



end
