function write_data_2_nc_technical(data, outfile, config)
% this function writes joyrad94 data into netcdf4

%################## Create a netCDF file.

ncid = netcdf.create(outfile,'NETCDF4'); 


%% ################## General settings for the fill values of the data

str_fill_value = 'NaN';
fill_value     = NaN  ;

%% ################# Define dimensions

did_time   = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
did_range  = netcdf.defDim(ncid,'range',data.n_levels);
did_vel    = netcdf.defDim(ncid,'velocity',max(data.DoppLen));
did_no_seq = netcdf.defDim(ncid,'number.chirp.sequences',data.no_chirp_seq);
did_scalar = netcdf.defDim(ncid,'scalar',1);

%% ######################## add global attributes
glob = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,glob,'FillValue',str_fill_value);
netcdf.putAtt(ncid,glob,'program_name',data.progname);
if data.DualPol == 0
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


% variables for dimensions
id_scal = netcdf.defVar(ncid,'scalar','nc_int',did_scalar);
netcdf.putAtt(ncid,id_scal,'comment','Variable for scalar-dimension');

id_no_seq = netcdf.defVar(ncid,'number.chirp.sequences','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_no_seq,'comment','Variable for number.chirp.sequences-dimension');

id_veldim = netcdf.defVar(ncid,'velocity','nc_int',did_vel);
netcdf.putAtt(ncid,id_veldim,'comment','Variable for velocity-dimension');



%% ################ get variable ids and add attributes

%%%%%%%%%% scalar variables

id_lat = netcdf.defVar(ncid,'lat','nc_float',did_scalar);
netcdf.putAtt(ncid,id_lat,'GEOMS_name','LATITUDE.INSTRUMENT');
netcdf.putAtt(ncid,id_lat,'long_name','latitude');
netcdf.putAtt(ncid,id_lat,'rpg_name','Lat');
netcdf.putAtt(ncid,id_lat,'standard_name','Latitude of instrument location (measurement site)');
netcdf.putAtt(ncid,id_lat,'short_name','lat');
netcdf.putAtt(ncid,id_lat,'units','degrees_north');
netcdf.putAtt(ncid,id_lat,'comment',['Latitude in degrees north [-90,90]']);

id_lon = netcdf.defVar(ncid,'lon','nc_float',did_scalar);
netcdf.putAtt(ncid,id_lon,'GEOMS_name','LONGITUDE.INSTRUMENT');
netcdf.putAtt(ncid,id_lon,'long_name','longitude');
netcdf.putAtt(ncid,id_lon,'rpg_name','Lon');
netcdf.putAtt(ncid,id_lon,'standard_name','Longitude of instrument location (measurement site)');
netcdf.putAtt(ncid,id_lon,'short_name','lon');
netcdf.putAtt(ncid,id_lon,'units','degrees_east');
netcdf.putAtt(ncid,id_lon,'comment',['Longitude in degrees east [-180,180]']);

id_MSL = netcdf.defVar(ncid,'zsl','nc_float',did_scalar);
netcdf.putAtt(ncid,id_MSL,'GEOMS_name','ALTITUDE.INSTRUMENT');
netcdf.putAtt(ncid,id_MSL,'long_name','altitude');
netcdf.putAtt(ncid,id_MSL,'rpg_name','MSL');
netcdf.putAtt(ncid,id_MSL,'standard_name','Altitude above mean sea level the instrument is located');
netcdf.putAtt(ncid,id_MSL,'short_name','zsl');
netcdf.putAtt(ncid,id_MSL,'units','m');
netcdf.putAtt(ncid,id_MSL,'comment',['Information in rpg_name (MSL) was added '...
                                     'during the data processing']);
                                 
id_freq = netcdf.defVar(ncid,'freq_sb','nc_float',did_scalar);
netcdf.putAtt(ncid,id_freq,'GEOMS_name','FREQUENCY');
netcdf.putAtt(ncid,id_freq,'long_name','radiation_frequency');
netcdf.putAtt(ncid,id_freq,'rpg_name','freq');
netcdf.putAtt(ncid,id_freq,'standard_name','Central transmission frequency of the electromagnetic wave emitted by radar');
netcdf.putAtt(ncid,id_freq,'short_name','freq_sb');
netcdf.putAtt(ncid,id_freq,'units','s-1');
netcdf.putAtt(ncid,id_freq,'comment',['Frequency can be converted into '...
                                      'WAVELENGHT = c/FREQUENCY; where c'...
                                      'is the speed of light']);

id_wl = netcdf.defVar(ncid,'wl','nc_float',did_scalar);
netcdf.putAtt(ncid,id_wl,'GEOMS_name','WAVELENGTH');
netcdf.putAtt(ncid,id_wl,'long_name','radiation_wavelength');
netcdf.putAtt(ncid,id_wl,'standard_name','Wavelength of the transmitted electromagnetic wave emitted by radar');
netcdf.putAtt(ncid,id_wl,'short_name','wl');
netcdf.putAtt(ncid,id_wl,'units','m');
netcdf.putAtt(ncid,id_wl,'comment',['Information in about operating radar '...
                                     'wavelengh was calculated based on the '...
                                     'frequency information.']);

id_HPBW = netcdf.defVar(ncid,'hpbw','nc_float',did_scalar);
netcdf.putAtt(ncid,id_HPBW,'GEOMS_name','ANTENNA_BEAM_WIDTH');
netcdf.putAtt(ncid,id_HPBW,'long_name','antenna_beam_width');
netcdf.putAtt(ncid,id_HPBW,'rpg_name','HPBW');
netcdf.putAtt(ncid,id_HPBW,'standard_name','Antenna half power beam width');
netcdf.putAtt(ncid,id_HPBW,'short_name','hpbw');
netcdf.putAtt(ncid,id_HPBW,'units','degrees');
netcdf.putAtt(ncid,id_HPBW,'comment',['Half power beam width is the angle between '...
                                      'the half-power (-3 dB) points of the main '...
                                      'lobe, when referenced to the peak effective '...
                                      'radiated power of the main lobe']);


id_AntiAlias = netcdf.defVar(ncid,'alias_flag','nc_byte',did_scalar);
netcdf.putAtt(ncid,id_AntiAlias,'GEOMS_name','processing.flag.alias_flag');
netcdf.putAtt(ncid,id_AntiAlias,'long_name','Quality flag for dealiasing');
netcdf.putAtt(ncid,id_AntiAlias,'comment',['The flag index shows: '...
                                           '0 = no dealiasing applied, '... 
                                           '1 = dealiasing by RPG, '...
                                           '2 = dealiasing by the applied code (see DATA_SOURCE)']);

%--- from here on the variables are not GEOMS formate conforme

id_filecode = netcdf.defVar(ncid,'filecode','nc_int',did_scalar);
netcdf.putAtt(ncid,id_filecode,'standard_name','filecode/version number');

id_samples = netcdf.defVar(ncid,'samples','nc_int',did_scalar);
netcdf.putAtt(ncid,id_samples,'standard_name','No. samples');

id_levels = netcdf.defVar(ncid,'n_levels','nc_int',did_scalar);
netcdf.putAtt(ncid,id_levels,'standard_name','Number of range levels');

id_no_chirp_seq = netcdf.defVar(ncid,'no_chirp_seq','nc_int',did_scalar);
netcdf.putAtt(ncid,id_no_chirp_seq,'standard_name','Number of chirp sequences');
netcdf.putAtt(ncid,id_no_chirp_seq,'comment',...
    'The radar can be programmed to run different resolution modes for different layers spanning several range gates');

id_CalInt = netcdf.defVar(ncid,'CalInt','nc_int',did_scalar);
netcdf.putAtt(ncid,id_CalInt,'standard_name','Sample interval for automated zero calibrations');
netcdf.putAtt(ncid,id_CalInt,'units','no. samples');

id_AntSep = netcdf.defVar(ncid,'AntSep','nc_float',did_scalar);
netcdf.putAtt(ncid,id_AntSep,'standard_name','Separation of antenna axis');
netcdf.putAtt(ncid,id_AntSep,'units','m');

id_AntDia = netcdf.defVar(ncid,'AntDia','nc_float',did_scalar);
netcdf.putAtt(ncid,id_AntDia,'standard_name','Antenna diameter');
netcdf.putAtt(ncid,id_AntDia,'units','m');

id_AntG = netcdf.defVar(ncid,'AntG','nc_float',did_scalar);
netcdf.putAtt(ncid,id_AntG,'standard_name','Antenna gain');
netcdf.putAtt(ncid,id_AntG,'units','[]');

id_SampDur = netcdf.defVar(ncid,'SampDur','nc_float',did_scalar);
netcdf.putAtt(ncid,id_SampDur,'standard_name','Full sample duration');
netcdf.putAtt(ncid,id_SampDur,'units','seconds');
netcdf.putAtt(ncid,id_SampDur,'comment','This is not the integration time (see SeqIntTime), but the sum of the chirp sequences integration times.')

id_C = netcdf.defVar(ncid,'C','nc_float',did_scalar);
netcdf.putAtt(ncid,id_C,'standard_name','Radar constant (see manual eq. 2.1.5)');

id_pol = netcdf.defVar(ncid,'DualPol','nc_byte',did_scalar);
netcdf.putAtt(ncid,id_pol,'standard_name','Polarisation. DualPol = 2: Full polarimetry. DualPol = 1: Polarimetric receiver but only vertical polarization transmitted. DualPol = 0: Only vertical polarization available.');

if isfield(data, 'CompEna')
    id_compress = netcdf.defVar(ncid,'CompEna','nc_byte',did_scalar);
    netcdf.putAtt(ncid,id_compress,'standard_name','If CompEna = 1, then compression was enabled in RPG software. Noise level from spectra was removed.');
end

id_cal_mom = netcdf.defVar(ncid,'cal_mom','nc_byte',did_scalar);
netcdf.putAtt(ncid,id_cal_mom,'standard_name','Integer indicating how moments were calculated.');
netcdf.putAtt(ncid,id_cal_mom,'comment',...
    'Set in config file. 2 = moments were calculated from spectra. 3 = moments are taken from RPG processing sofware (lv1).');



%%%%%%% range variables

id_range = netcdf.defVar(ncid,'range','nc_float',did_range);
netcdf.putAtt(ncid,id_range,'GEOMS_name','RANGE');
netcdf.putAtt(ncid,id_range,'long_name','range');
netcdf.putAtt(ncid,id_range,'rpg_name','range');
netcdf.putAtt(ncid,id_range,'standard_name','Range gates of the radar measurements');
netcdf.putAtt(ncid,id_range,'short_name','range');
netcdf.putAtt(ncid,id_range,'units','m');
netcdf.putAtt(ncid,id_range,'valid_range', [nanmin(data.range(:)), nanmax(data.range(:))]);
netcdf.putAtt(ncid,id_range,'comment',['Range from the radar antenna to the '...
                                       'center of the radar range gate']);
                                   

%--- from here on the variables are not GEOMS formate conforme

id_Fr = netcdf.defVar(ncid,'Fr','nc_int',did_range);
netcdf.putAtt(ncid,id_Fr,'standard_name','Range factors, see manual eq. (2.5.6)');

id_range_offsets = netcdf.defVar(ncid,'range_offsets','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_range_offsets,'GEOMS_name','radar.operation.parameter.range_index_chirp');
netcdf.putAtt(ncid,id_range_offsets,'long_name','range_index_of_chirp_sequence_change');
netcdf.putAtt(ncid,id_range_offsets,'rpg_name','range_offsets');
netcdf.putAtt(ncid,id_range_offsets,'standard_name','Index indicating the range gate where the new/next Chirp sequence starts');
netcdf.putAtt(ncid,id_range_offsets,'units','');
netcdf.putAtt(ncid,id_range_offsets,'comment',['The command range(range_offsets) '...
                                               'will give you the range where a '...
                                               'new chirp sequence starts. '...
                                               'range_offsets counts from '...
                                               '1 to n_levels.']);

                                           
%%%%%%%% chirp_seq_dependent variables

id_DoppMax = netcdf.defVar(ncid,'nqv','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_DoppMax,'GEOMS_name','NYQUIST.VELOCITY');
netcdf.putAtt(ncid,id_DoppMax,'long_name','nyquist_velocity');
netcdf.putAtt(ncid,id_DoppMax,'rpg_name','DoppMax');
netcdf.putAtt(ncid,id_DoppMax,'standard_name',['Max. unambigious Doppler velocity for each chirp sequence; '...
                                          'Nyquist velocity per chirp sequence']);
netcdf.putAtt(ncid,id_DoppMax,'short_name','nqv');
netcdf.putAtt(ncid,id_DoppMax,'units','m s-1');
netcdf.putAtt(ncid,id_DoppMax,'comment',['The Nyquist velocity is the maximum velocity '...
                                         'that can be correctly displayed by a Doppler '...
                                         'radar, and this is dependent on the wavelength '...
                                         'and frequency of eleoctromagnetic wave '...
                                         'emitted by the radar.']);

id_SeqAvg = netcdf.defVar(ncid,'seq_avg','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_SeqAvg,'GEOMS_name','radar.operation.parameter.avg_chirps_per_chirp');
netcdf.putAtt(ncid,id_SeqAvg,'standard_name','Number of averaged chirps in each chirp sequence');

id_SeqIntTime = netcdf.defVar(ncid,'seq_int_time','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_SeqIntTime,'GEOMS_name','radar.operation.parameter.int_time_chirp');
netcdf.putAtt(ncid,id_SeqIntTime,'standard_name','Integration time of each chirp sequence');
netcdf.putAtt(ncid,id_SeqIntTime,'units','seconds');

id_DoppLen = netcdf.defVar(ncid,'dopp_len','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_DoppLen,'GEOMS_name','radar.operation.parameter.spec_samples_per_chirp');
netcdf.putAtt(ncid,id_DoppLen,'standard_name','Number of samples in Doppler spectra of each chirp sequence');
netcdf.putAtt(ncid,id_DoppLen,'comment','Needed to calculate the Doppler resolution: DoppRes = 2*nqv/dopp_len');

id_nAvg = netcdf.defVar(ncid,'n_avg','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_nAvg,'GEOMS_name','radar.operation.parameter.no_spec_avg');
netcdf.putAtt(ncid,id_nAvg,'standard_name','Number of spectra averaged for each chirp');
netcdf.putAtt(ncid,id_nAvg,'comment','n_avg = seq_avg/dopp_len')


%--- from here on the variables are not GEOMS formate conforme

id_ChrpIntTime = netcdf.defVar(ncid,'chirp_seq_int_time','nc_float',did_no_seq); % added by RG 4.12.2020
netcdf.putAtt(ncid,id_ChrpIntTime,'standard_name',['Integration time of each chirp sequence, ' ...
                                                  'calculated from chirp repetition frequency ' ...
                                                  'and number of repetitions']);
netcdf.putAtt(ncid,id_ChrpIntTime,'units','seconds');
netcdf.putAtt(ncid,id_ChrpIntTime,'comment',['This variable corresponds to the chirp ' ...
                                             'integration time set in the chirp table ' ...
                                             'during mdf definition']);
                                         
%%%%%%%% time dependend variables

id_time = netcdf.defVar(ncid,'time','nc_uint',did_time);
netcdf.putAtt(ncid,id_time,'GEOMS_name','DATETIME');
netcdf.putAtt(ncid,id_time,'long_name','time');
netcdf.putAtt(ncid,id_time,'rpg_name','time');
netcdf.putAtt(ncid,id_time,'standard_name','time stamp of the radar measurements');
netcdf.putAtt(ncid,id_time,'short_name','time');
netcdf.putAtt(ncid,id_time,'units','seconds since 2001-1-1 0:0:0');
netcdf.putAtt(ncid,id_time,'valid_range', [nanmin(data.time(:)), nanmax(data.time(:))]);
netcdf.putAtt(ncid,id_time,'comment','Time in UTC');
if isfield(data, 'totsampchangelabel' )
    netcdf.putAtt(ncid,id_time, 'quality_flag', 'Dublicate time stamps found in lv0-file, the first occurrence of the dublicate time is removed')
end



id_status = netcdf.defVar(ncid,'blower_status','nc_byte',did_time);
netcdf.putAtt(ncid,id_status,'GEOMS_name','BLOWER.STATUS');
netcdf.putAtt(ncid,id_status,'long_name','platform_blower_status');
netcdf.putAtt(ncid,id_status,'rpg_name','status');
netcdf.putAtt(ncid,id_status,'standard_name','Status of the blower which keeps the radar antennas dry');
netcdf.putAtt(ncid,id_status,'short_name','blowe_status');
netcdf.putAtt(ncid,id_status,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_status,'comment',['The quality flag shows the following '...
                                        'status: 0/1 = heater on/off; '...
                                        '0/10 = blower on/off. The parameter '...
                                        'is recorded to check and monitor '...
                                        'the radar performance']);

id_TransPow = netcdf.defVar(ncid,'p_trans','nc_float',did_time);
netcdf.putAtt(ncid,id_TransPow,'GEOMS_name','TRANSMITTER.POWER');
netcdf.putAtt(ncid,id_TransPow,'long_name','platform_transmitter_power');
netcdf.putAtt(ncid,id_TransPow,'rpg_name','TransPow');
netcdf.putAtt(ncid,id_TransPow,'standard_name','Transmitted power of the radar instrument');
netcdf.putAtt(ncid,id_TransPow,'short_name','p_trans');
netcdf.putAtt(ncid,id_TransPow,'units','W');
netcdf.putAtt(ncid,id_TransPow,'valid_range',[nanmin(data.TransPow(:)), nanmax(data.TransPow(:))]);
netcdf.putAtt(ncid,id_TransPow,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_TransPow,'comment',['The transmitted power of the radar instrument '...
                                          'is part of the house keeping data and '...
                                          'used to monitor its the performance.']);
                                 
id_T_trans = netcdf.defVar(ncid,'t_trans','nc_float',did_time);
netcdf.putAtt(ncid,id_T_trans,'GEOMS_name','TRANSMITTER.TEMPERATURE');
netcdf.putAtt(ncid,id_T_trans,'long_name','platform_transmitter_temperature');
netcdf.putAtt(ncid,id_T_trans,'rpg_name','T_trans');
netcdf.putAtt(ncid,id_T_trans,'standard_name','Transmitter temperature of the radar instrument');
netcdf.putAtt(ncid,id_T_trans,'short_name','t_trans');
netcdf.putAtt(ncid,id_T_trans,'units','K');
netcdf.putAtt(ncid,id_T_trans,'valid_range',[nanmin(data.T_trans(:)), nanmax(data.T_trans(:))]);
netcdf.putAtt(ncid,id_T_trans,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_T_trans,'comment',['The transmitter temperature of the radar '...
                                         'instrument is part of the house keeping data '...
                                         'and used to monitor its the performance. '...
                                         'For a stable performance 309K < t_trans < 313K.']);
                                 
id_T_rec = netcdf.defVar(ncid,'t_rec','nc_float',did_time);
netcdf.putAtt(ncid,id_T_rec,'GEOMS_name','RECEIVER.TEMPERATURE');
netcdf.putAtt(ncid,id_T_rec,'long_name','platform_receiver_temperature');
netcdf.putAtt(ncid,id_T_rec,'rpg_name','T_rec');
netcdf.putAtt(ncid,id_T_rec,'standard_name','Receiver temperature of the radar instrument');
netcdf.putAtt(ncid,id_T_rec,'short_name','t_rec');
netcdf.putAtt(ncid,id_T_rec,'units','K');
netcdf.putAtt(ncid,id_T_rec,'valid_range',[nanmin(data.T_rec(:)), nanmax(data.T_rec(:))]);
netcdf.putAtt(ncid,id_T_rec,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_T_rec,'comment',['The receiver temperature of the radar '...
                                       'instrument is part of the house keeping data '...
                                       'and used to monitor its the performance. '...
                                       'For a stable performance 306K < t_rec < 312K.']);
                                 
id_T_pc = netcdf.defVar(ncid,'t_pc','nc_float',did_time);
netcdf.putAtt(ncid,id_T_pc,'GEOMS_name','COMPUTER.TEMPERATURE');
netcdf.putAtt(ncid,id_T_pc,'long_name','platform_computer_temperature');
netcdf.putAtt(ncid,id_T_pc,'rpg_name','T_pc');
netcdf.putAtt(ncid,id_T_pc,'standard_name','Temperature of the radar controlling PC');
netcdf.putAtt(ncid,id_T_pc,'short_name','t_pc');
netcdf.putAtt(ncid,id_T_pc,'units','K');
netcdf.putAtt(ncid,id_T_pc,'valid_range',[nanmin(data.T_pc(:)), nanmax(data.T_pc(:))]);
netcdf.putAtt(ncid,id_T_pc,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_T_pc,'comment',['The radar controlling PC temperature '...
                                     'is part of the house keeping data and '...
                                     'used to monitor its the performance. '...
                                     'For a stable performance t_pc < 323K.']);

id_QF = netcdf.defVar(ncid,'rad_qf','nc_byte',did_time);
netcdf.putAtt(ncid,id_QF,'GEOMS_name','FLAG.MEASUREMENT.QUALITY');
netcdf.putAtt(ncid,id_QF,'long_name','quality_falg_platform_measurement');
netcdf.putAtt(ncid,id_QF,'rpg_name','QF');
netcdf.putAtt(ncid,id_QF,'standard_name',['Quality flag related to the radar operation '... 
                                      'given by radar software from RPG']);
netcdf.putAtt(ncid,id_QF,'short_name','rad_qf');                                  
netcdf.putAtt(ncid,id_QF,'comment', ['To get the bit entries, one has to'...
                                     'convert the integer into a 4 bit binary. '...
                                     'bit4 = ADC saturation, bit3 = spectral '...
                                     'width too high, bit2 = no transmitter '...
                                     'power leveling. Note that in the above '...
                                     'convention holds: bit1 = 2^3, '...
                                     'bit2 = 2^2, bit3 = 2^1, bit4 = 2^0'])

id_sampleTms = netcdf.defVar(ncid,'sample_tms','nc_int',did_time);
netcdf.putAtt(ncid,id_sampleTms,'GEOMS_name','DATETIME.MILLISECONDS');
netcdf.putAtt(ncid,id_sampleTms,'standard_name','Milliseconds of sample');
netcdf.putAtt(ncid,id_sampleTms,'units','mu s');
netcdf.putAtt(ncid,id_sampleTms,'comment','To get the correct time the variable sample_tms must be added: time = time + sample_tms.');

id_incel = netcdf.defVar(ncid,'ins_elevation','nc_float',did_time); 
netcdf.putAtt(ncid,id_incel,'GEOMS_name','PLATFORM.ANGLE.ELEVATION');
netcdf.putAtt(ncid,id_incel,'long_name','platform_roll_angle');
netcdf.putAtt(ncid,id_incel,'standard_name','Inclination angle of the instrument; elevation');
netcdf.putAtt(ncid,id_incel,'units','deg');
                           
id_incea = netcdf.defVar(ncid,'ins_azimuth','nc_float',did_time);
netcdf.putAtt(ncid,id_incea,'GEOMS_name','PLATFORM.ANGLE.AZIMUTH');
netcdf.putAtt(ncid,id_incea,'long_name','angle_of_rotation_from_solar_azimuth_to_platform_azimuth');
netcdf.putAtt(ncid,id_incea,'standard_name','Inclination angle of the instrument; azimuth');
netcdf.putAtt(ncid,id_incea,'units','deg');

id_ele = netcdf.defVar(ncid,'ele','nc_float',did_time);
netcdf.putAtt(ncid,id_ele,'GEOMS_name','operation.parameter.system_elevation');
netcdf.putAtt(ncid,id_ele,'standard_name','Sensor elevation angle');
netcdf.putAtt(ncid,id_ele,'units','degrees');

id_az = netcdf.defVar(ncid,'azi','nc_float',did_time);
netcdf.putAtt(ncid,id_az,'GEOMS_name','operation.parameter.system_azimuth');
netcdf.putAtt(ncid,id_az,'standard_name','Sensor azimuth angle');
netcdf.putAtt(ncid,id_az,'units','degrees');

id_vol = netcdf.defVar(ncid,'vol','nc_float',did_time);
netcdf.putAtt(ncid,id_vol,'GEOMS_name','operation.parameter.vol');
netcdf.putAtt(ncid,id_vol,'standard_name','Voltage direct detection channel');
netcdf.putAtt(ncid,id_vol,'units','V');

id_powIF = netcdf.defVar(ncid,'pow_if','nc_float',did_time);
netcdf.putAtt(ncid,id_powIF,'GEOMS_name','operation.parameter.pow_if');
netcdf.putAtt(ncid,id_powIF,'standard_name','IF power at ADC');
netcdf.putAtt(ncid,id_powIF,'units','microWatts');
                                 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% multi-D variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

id_Aliasmask = netcdf.defVar(ncid,'alias_mask','nc_byte',[did_range,did_time]);
netcdf.putAtt(ncid,id_Aliasmask,'GEOMS_name','MASK.PROCESSING.ALIAISING');
netcdf.putAtt(ncid,id_Aliasmask,'long_name','alaising_mask');
netcdf.putAtt(ncid,id_Aliasmask,'standard_name','Mask array indicating in which bin aliasing occurs');
netcdf.putAtt(ncid,id_Aliasmask,'comment',['AliasMask indicates the bins where aliasing is detected. '...
                                           '0 = no aliasing; 1 = aliasing occurs. Only applicable if '...
                                           'variable alias_flag is 2.']);
         
id_AliasStatus = netcdf.defVar(ncid,'AliasStatus','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_AliasStatus,'standard_name','Flags indicating quality of dealiasing.');
netcdf.putAtt(ncid,id_AliasStatus,'comment',['Must be converted into four bit binary string. If 0, i.e. dec2bin(AliasStatus) = 0000 then dealiasing was successful. ',...
    'If 2^0 bit is 1, there was no initial guess velocity for this bin. ',...
    'If 2^1 bit is 1, the spectrum where the main peak was assumed to occur could not be concatenated properly since chirp sequence boundaries were reached. ',...
    'If 2^2 bit is 1, there is still significant signal close to the nyquist limit. dealiasing might not have been performed properly. ',...
    'If 2^3 bit is 1, the mean of the mean Doppler velocity in differs by more than 5 m/s from the neighbouring bins mean value. ',...
    'The larger the decimal number, the less reliable the data. ',...
    'This flag is only available if dealiasing was applied in process_joyrad94_data.m. For RPG dealising there is no quality indicator.']);



id_vel = netcdf.defVar(ncid,'chirp_vel','nc_float',[did_vel,did_no_seq]);
netcdf.putAtt(ncid,id_vel,'standard_name','Velocity arrays for each chirp sequence');
netcdf.putAtt(ncid,id_vel,'units','m/s');
netcdf.putAtt(ncid,id_vel,'comment',...
    ['If spectra has been dealiased, use MinVel to get the correct velocity array for each spectra.' ...
    'The velocity array is asymmetric, i.e. the absolute values of maximum and minumum velocities are not equal. ',...
    'Since the spectrum at -v_nyquist and +v_nyquist is the same, the entry at +v_nyquist was cut.']);

id_MinVel = netcdf.defVar(ncid,'MinVel','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_MinVel,'standard_name','Minimum velocity in velocity array');
netcdf.putAtt(ncid,id_MinVel,'comment',['If dealising was applied the array indicates the velocity of the first bin in spec. ',...
    'The correct velocity array can be obtained by: vel_true = chirp_vel + Minvel(i,j) - chirp_vel(1)']);
netcdf.putAtt(ncid,id_MinVel,'units','m/s');

id_MinVel_Correction = netcdf.defVar(ncid,'MinVel_Correction','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_MinVel_Correction,'standard_name','Minimum velocity correction');
netcdf.putAtt(ncid,id_MinVel_Correction,'comment',['If dealising was applied the array indicates by how much MinVel '...
    'was corrected by the final quality check (output from dealias_spectra_vm_column_quality_check.m). '...
    'Subtracting this value from MinVel provides the offset before the final quality check (output from dealias_spectra.m).']);
netcdf.putAtt(ncid,id_MinVel_Correction,'units','m/s');

id_mask = netcdf.defVar(ncid,'mask','nc_byte',[did_range,did_time]);
netcdf.putAtt(ncid,id_mask,'standard_name','Mask array of occupied range cells: 0 = not occupied; 1 = occupied');
netcdf.putAtt(ncid,id_mask,'comment',['Mask(time,range) is set to 0 when signal does not exceed threshold set in the software. ',...
    'The threshold is a multiple of the noise standard deviation that is set in the measurement definition file and defines ',...
    'the linear sensitivity limit (see variable SLv). Doppler spectrum is only stored when set to 1. ',...
    'Note that data produced with software version 1 (i.e. before 7 Feb 2017) used a different threshold, ',...
    'i.e. the standard deviation for each chrip sequence. ',...
    'Moreover, the latter was calculated wrongly by the software version 1, which lead to the rejection of actual valid samples.']);


id_PNv = netcdf.defVar(ncid,'PNv','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_PNv,'standard_name','Total IF power in vertical polarization measured at ADC input');
netcdf.putAtt(ncid,id_PNv,'units','W');

id_VNoisePow_peak = netcdf.defVar(ncid,'v_noise_pow_peak','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_VNoisePow_peak,'GEOMS_name','radar.doppler.spectrum_noise.power.peak_vv');
netcdf.putAtt(ncid,id_VNoisePow_peak,'standard_name','Bin Doppler spectrum peak noise power in vertical polarization');
netcdf.putAtt(ncid,id_VNoisePow_peak,'units','mm6 m-3');
netcdf.putAtt(ncid,id_VNoisePow_peak,'valid_range',[nanmin(data.VNoisePow_peak(:)), nanmax(data.VNoisePow_peak(:))]);
netcdf.putAtt(ncid,id_VNoisePow_peak,'fill_value','NaNf');
netcdf.putAtt(ncid,id_VNoisePow_peak,'comment','If cal_mom == 3, then bin noise power is only calculated if CompEna == 1.');


if data.DualPol > 0
    
    id_PNh = netcdf.defVar(ncid,'if_pow_hh','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_PNh,'standard_name','Total IF power in horizontal polarization measured at ADC input');
    netcdf.putAtt(ncid,id_PNh,'units','W');
    
end % if data.DualPol > 0


if data.CompEna > 0 && data.DualPol > 0
         id_HNoisePow_peak = netcdf.defVar(ncid,'h_noise_pow_peak','nc_float',[did_range,did_time]);
         netcdf.putAtt(ncid,id_HNoisePow_peak,'GEOMS_name','radar.doppler.spectrum_noise.power.peak_hh');                  
         netcdf.putAtt(ncid,id_HNoisePow_peak,'standard_name','Bin Doppler spectrum peak noise power in horizontal polarization');
         netcdf.putAtt(ncid,id_HNoisePow_peak,'units','mm6 mm-3');
         netcdf.putAtt(ncid,id_HNoisePow_peak,'comment','If cal_mom == 3, then bin noise power is only calculated if CompEna == 1.');
end % data.CompEna > 0



%% ###################### initialize compression

netcdf.defVarDeflate(ncid,id_vol,true,true,9);
netcdf.defVarDeflate(ncid,id_powIF,true,true,9);
netcdf.defVarDeflate(ncid,id_status,true,true,9);
netcdf.defVarDeflate(ncid,id_TransPow,true,true,9);
netcdf.defVarDeflate(ncid,id_T_trans,true,true,9);
netcdf.defVarDeflate(ncid,id_T_rec,true,true,9);
netcdf.defVarDeflate(ncid,id_T_pc,true,true,9);
netcdf.defVarDeflate(ncid,id_QF,true,true,9);
netcdf.defVarDeflate(ncid,id_incel,true,true,9);
netcdf.defVarDeflate(ncid,id_incea,true,true,9);
netcdf.defVarDeflate(ncid,id_ele,true,true,9);
netcdf.defVarDeflate(ncid,id_az,true,true,9);
netcdf.defVarDeflate(ncid,id_mask,true,true,9);
netcdf.defVarDeflate(ncid,id_Aliasmask,true,true,9);
netcdf.defVarDeflate(ncid,id_AliasStatus,true,true,9);
netcdf.defVarDeflate(ncid,id_MinVel,true,true,9);
netcdf.defVarDeflate(ncid,id_PNv,true,true,9);
netcdf.defVarDeflate(ncid,id_VNoisePow_peak,true,true,9);
netcdf.defVarDeflate(ncid,id_MinVel_Correction,true,true,9);

if data.CompEna > 0 && data.DualPol > 0
    netcdf.defVarDeflate(ncid,id_HNoisePow_peak,true,true,9);
end    

% variables for dimensions
netcdf.defVarDeflate(ncid,id_scal,true,true,9);
netcdf.defVarDeflate(ncid,id_no_seq,true,true,9);
netcdf.defVarDeflate(ncid,id_veldim,true,true,9);


netcdf.endDef(ncid);


%% ####################### put variables into file

% variables for dimensions
netcdf.putVar(ncid,id_scal,0,1);
netcdf.putVar(ncid,id_no_seq,0,data.no_chirp_seq,1:data.no_chirp_seq);
netcdf.putVar(ncid,id_veldim,0,max(data.DoppLen), 1:max(data.DoppLen));


% scalars
netcdf.putVar(ncid,id_lat,0,data.Lat);
netcdf.putVar(ncid,id_lon,0,data.Lon);
netcdf.putVar(ncid,id_MSL,0,config.MSL);
netcdf.putVar(ncid,id_freq,0,data.freq);
netcdf.putVar(ncid, id_wl, 0, data.freq * 1e9);
netcdf.putVar(ncid,id_HPBW,0,data.HPBW);
netcdf.putVar(ncid,id_AntiAlias,0,data.AntiAlias);
netcdf.putVar(ncid,id_filecode,0,data.filecode);
netcdf.putVar(ncid,id_samples,0,data.totsamp);
netcdf.putVar(ncid,id_levels,0,data.n_levels);
netcdf.putVar(ncid,id_no_chirp_seq,0,data.no_chirp_seq);
netcdf.putVar(ncid,id_CalInt,0,data.CalInt);
netcdf.putVar(ncid,id_AntSep,0,data.AntSep);
netcdf.putVar(ncid,id_AntDia,0,data.AntDia);
netcdf.putVar(ncid,id_AntG,0,data.AntG);
netcdf.putVar(ncid,id_SampDur,0,data.SampDur);
netcdf.putVar(ncid,id_C,0,data.C);
netcdf.putVar(ncid,id_pol,0,data.DualPol);
if isfield(data, 'CompEna')
    netcdf.putVar(ncid,id_compress,0,data.CompEna);
end
netcdf.putVar(ncid,id_cal_mom,0,data.cal_mom);


% range dependent
netcdf.putVar(ncid,id_range,0,data.n_levels,data.range);
netcdf.putVar(ncid,id_Fr,0,data.n_levels,data.Fr);

% chrip seq dependent variables
netcdf.putVar(ncid,id_DoppMax,0,data.no_chirp_seq,data.DoppMax);
netcdf.putVar(ncid,id_SeqAvg,0,data.no_chirp_seq,data.SeqAvg);
netcdf.putVar(ncid,id_SeqIntTime,0,data.no_chirp_seq,data.SeqIntTime);
netcdf.putVar(ncid,id_DoppLen,0,data.no_chirp_seq,data.DoppLen);
netcdf.putVar(ncid,id_nAvg,0,data.no_chirp_seq,data.nAvg);
netcdf.putVar(ncid,id_ChrpIntTime,0,data.no_chirp_seq,data.ChirpIntTime); % added by RG 4.12.2020

% time dependent variables
netcdf.putVar(ncid,id_time,0,data.totsamp,data.time);
netcdf.putVar(ncid,id_status,0,data.totsamp,data.status);
netcdf.putVar(ncid,id_TransPow,0,data.totsamp,data.TransPow);
netcdf.putVar(ncid,id_T_trans,0,data.totsamp,data.T_trans);
netcdf.putVar(ncid,id_T_rec,0,data.totsamp,data.T_rec);
netcdf.putVar(ncid,id_T_pc,0,data.totsamp,data.T_pc);
netcdf.putVar(ncid,id_QF,0,data.totsamp,data.QF);
netcdf.putVar(ncid,id_sampleTms,0,data.totsamp,data.sampleTms);
netcdf.putVar(ncid,id_incel,0,data.totsamp,data.reserved(:,2));
netcdf.putVar(ncid,id_incea,0,data.totsamp,data.reserved(:,3));
netcdf.putVar(ncid,id_ele,0,data.totsamp,data.ele);
netcdf.putVar(ncid,id_az,0,data.totsamp,data.az);
netcdf.putVar(ncid,id_vol,0,data.totsamp,data.u);
netcdf.putVar(ncid,id_powIF,0,data.totsamp,data.powIF);

% multidimensional variables
netcdf.putVar(ncid,id_Aliasmask,[0,0],[data.n_levels,data.totsamp],data.Aliasmask');
netcdf.putVar(ncid,id_AliasStatus,[0,0],[data.n_levels,data.totsamp],data.AliasStatus');
netcdf.putVar(ncid,id_vel,[0,0],[max(data.DoppLen),data.no_chirp_seq],data.velocity');
netcdf.putVar(ncid,id_MinVel,[0,0],[data.n_levels,data.totsamp],data.MinVel');
netcdf.putVar(ncid,id_MinVel_Correction,[0,0],[data.n_levels,data.totsamp],data.MinVel_Correction');
netcdf.putVar(ncid,id_mask,[0,0],[data.n_levels,data.totsamp],data.mask');
netcdf.putVar(ncid,id_PNv,[0,0],[data.n_levels,data.totsamp],data.PNv');
netcdf.putVar(ncid,id_VNoisePow_peak,[0,0],[data.n_levels,data.totsamp],data.VNoisePow_peak');


if data.DualPol > 0
    netcdf.putVar(ncid,id_PNh,[0,0],[data.n_levels,data.totsamp],data.PNh');    
end


if data.CompEna > 0 && data.DualPol > 0
    netcdf.putVar(ncid,id_HNoisePow_peak,[0,0],[data.n_levels,data.totsamp],data.HNoisePow_peak');
end

netcdf.close(ncid);
