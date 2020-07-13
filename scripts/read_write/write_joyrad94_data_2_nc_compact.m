function write_joyrad94_data_2_nc_compact(data, outfile, config)

% this function writes joyrad94 data into netcdf4
% Changes of all the 

%% ################## Create a netCDF file.

ncid = netcdf.create(outfile,'NETCDF4'); 


%% ################## General settings for the fill values of the data

str_fill_value = '-999';
fill_value     = -999  ;

%% ################# Define dimensions

%did_time = netcdf.defDim(ncid,'time',data.totsamp);
%Changed by J.A. Bravo-Aranda
did_time = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
did_range = netcdf.defDim(ncid,'range',data.n_levels);
did_no_seq = netcdf.defDim(ncid,'number.chirp.sequences',data.no_chirp_seq);
did_scalar = netcdf.defDim(ncid,'scalar',1);


%% ######################## add global attributes
glob = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,glob,'fill_value',str_fill_value);
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
netcdf.putAtt(ncid,id_lon,'standard_name','Longitude of instrument locarion (measurement site)');
netcdf.putAtt(ncid,id_lon,'short_name','lon');
netcdf.putAtt(ncid,id_lon,'units','degrees_east');
netcdf.putAtt(ncid,id_lon,'comment',['Longitude in degrees east [-180,180]']);

id_MSL = netcdf.defVar(ncid,'zsl','nc_float',did_scalar);
netcdf.putAtt(ncid,id_MSL,'GEOMS_name','ALTITUDE.INSTRUMENT');
netcdf.putAtt(ncid,id_MSL,'long_name','altitude');
netcdf.putAtt(ncid,id_MSL,'rpg_name','MSL');
netcdf.putAtt(ncid,id_MSL,'standard_name','Altitud above mean sea level the instrument is located');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% range variables %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

id_range = netcdf.defVar(ncid,'range','nc_float',did_range);
netcdf.putAtt(ncid,id_range,'GEOMS_name','RANGE');
netcdf.putAtt(ncid,id_range,'long_name','range');
netcdf.putAtt(ncid,id_range,'rpg_name','range');
netcdf.putAtt(ncid,id_range,'standard_name','Range gates of the radar measuremnts');
netcdf.putAtt(ncid,id_range,'short_name','range');
netcdf.putAtt(ncid,id_range,'units','m');
netcdf.putAtt(ncid,id_range,'valid_range', [nanmin(data.range(:)), nanmax(data.range(:))]);
netcdf.putAtt(ncid,id_range,'comment',['Range from the radar antenna to the '...
                                       'center of the radar range gate']);


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

% id_AntiAlias = netcdf.defVar(ncid,'alias_flag','nc_byte',did_scalar);
% netcdf.putAtt(ncid,id_AntiAlias,'GEOMS_name','processing.flag.alias_flag');
% netcdf.putAtt(ncid,id_AntiAlias,'long_name','Quality flag for dealiasing');
% netcdf.putAtt(ncid,id_AntiAlias,'comment',['The falg index shows: '...
%                                            '0 = no dealiasing applied, '... 
%                                            '1 = dealiasing by RPG, '...
%                                            '2 = dealiasing by the applied code (see DATA_SOURCE)']);


%%%%%%%% time dependend variables

id_time = netcdf.defVar(ncid,'time','nc_uint',did_time);
netcdf.putAtt(ncid,id_time,'GEOMS_name','DATETIME');
netcdf.putAtt(ncid,id_time,'long_name','time');
netcdf.putAtt(ncid,id_time,'rpg_name','time');
netcdf.putAtt(ncid,id_time,'standard_name','time stap of the radar measurments');
netcdf.putAtt(ncid,id_time,'short_name','time');
netcdf.putAtt(ncid,id_time,'units','Time in sec since 2001.01.01. 00:00:00.');
netcdf.putAtt(ncid,id_time,'valid_range', [nanmin(data.time(:)), nanmax(data.time(:))]);
netcdf.putAtt(ncid,id_time,'comment','Time formate is MJD2K');
if isfield(data, 'totsampchangelabel' )
    netcdf.putAtt(ncid,id_time, 'quality_flag', 'Dublicate time stamps found in lv0-file, the first occurrence of the dublicate time is removed')
end

% id_sampleTms = netcdf.defVar(ncid,'sample_tms','nc_int',did_time);
% netcdf.putAtt(ncid,id_sampleTms,'GEOMS_name','DATETIME.milliseconds');
% netcdf.putAtt(ncid,id_sampleTms,'long_name','Milliseconds of sample');
% netcdf.putAtt(ncid,id_sampleTms,'units','mu s');
% netcdf.putAtt(ncid,id_sampleTms,'comment','To get the correct time the variable sample_tms must be added: time = time + sample_tms.');

id_RR = netcdf.defVar(ncid,'rr','nc_float',did_time);
netcdf.putAtt(ncid,id_RR,'GEOMS_name','RAIN.RATE.SURFACE');
netcdf.putAtt(ncid,id_RR,'long_name','rainfall_rate');
netcdf.putAtt(ncid,id_RR,'rpg_name','RR');
netcdf.putAtt(ncid,id_RR,'standard_name','Rain rate measured by the meteo-station attached to the radar');
netcdf.putAtt(ncid,id_RR,'short_name','rr');
netcdf.putAtt(ncid,id_RR,'units','mm h-1');
netcdf.putAtt(ncid,id_RR,'valid_range',[nanmin(data.RR(:)), nanmax(data.RR(:))]);
netcdf.putAtt(ncid,id_RR,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_RR,'RAIN.RATE.SURFACE_SOURCE','Vaisala weather station WXT520 or WXT530');

id_rh = netcdf.defVar(ncid,'rh','nc_float',did_time);
netcdf.putAtt(ncid,id_rh,'GEOMS_name','HUMIDITY.RELATIVE.SURFACE');
netcdf.putAtt(ncid,id_rh,'long_name','relative_humidiry');
netcdf.putAtt(ncid,id_rh,'rpg_name','rh');
netcdf.putAtt(ncid,id_rh,'standard_name','Relative humidity measured by the meteo-station attached to the radar');
netcdf.putAtt(ncid,id_rh,'short_name','rh');
netcdf.putAtt(ncid,id_rh,'units','%');
netcdf.putAtt(ncid,id_rh,'valid_range',[nanmin(data.rh(:)), nanmax(data.rh(:))]);
netcdf.putAtt(ncid,id_rh,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_rh,'HUMIDITY.RELATIVE.SURFACE_SOURCE','Vaisala weather station WXT520 or WXT530');

id_T_env = netcdf.defVar(ncid,'ta','nc_float',did_time);
netcdf.putAtt(ncid,id_T_env,'GEOMS_name','TEMPERATURE.SURFACE');
netcdf.putAtt(ncid,id_T_env,'long_name','surface_temperature');
netcdf.putAtt(ncid,id_T_env,'rpg_name','T_env');
netcdf.putAtt(ncid,id_T_env,'standard_name','Temperature of the environment measured by the meteo-station attached to the radar');
netcdf.putAtt(ncid,id_T_env,'short_name','ta');
netcdf.putAtt(ncid,id_T_env,'units','K');
netcdf.putAtt(ncid,id_T_env,'valid_range',[nanmin(data.T_env(:)), nanmax(data.T_env(:))]);
netcdf.putAtt(ncid,id_T_env,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_T_env,'TEMPERATURE.SURFACE_SOURCE','Vaisala weather station WXT520 or WXT530');
netcdf.putAtt(ncid,id_T_env,'comment',['Air temperature is the bulk temperature '...
                                       'of the air, not the surface (skin) temperature.']);

id_pres = netcdf.defVar(ncid,'pa','nc_float',did_time);
netcdf.putAtt(ncid,id_pres,'GEOMS_name','SURFACE.PRESSURE');
netcdf.putAtt(ncid,id_pres,'long_name','surface_air_pressure');
netcdf.putAtt(ncid,id_pres,'rpg_name','press');
netcdf.putAtt(ncid,id_pres,'standard_name','Surface air pressure measured by the meteo-station attached to the radar');
netcdf.putAtt(ncid,id_pres,'short_name','pa');
netcdf.putAtt(ncid,id_pres,'units','hPa');
netcdf.putAtt(ncid,id_pres,'valid_range',[nanmin(data.pres(:)), nanmax(data.pres(:))]);
netcdf.putAtt(ncid,id_pres,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_pres,'SURFACE.PRESSURE_SOURCE','Vaisala weather station WXT520 or WXT530');
netcdf.putAtt(ncid,id_pres,'comment',['The surface called "surface" means the '...
                                      'lower boundary of the atmosphere. Air '...
                                      'pressure is the force per unit area '...
                                      'which would be exerted when the moving '...
                                      'gas molecules of which the air is '...
                                      'composed strike a theoretical surface '...
                                      'of any orientation.']);

id_ff = netcdf.defVar(ncid,'wspeed','nc_float',did_time);
netcdf.putAtt(ncid,id_ff,'GEOMS_name','WIND.SPEED.SURFACE');
netcdf.putAtt(ncid,id_ff,'long_name','wind_speed');
netcdf.putAtt(ncid,id_ff,'rpg_name','ff');
netcdf.putAtt(ncid,id_ff,'standard_name','wind speed measured by the meteo-station attached to the radar');
netcdf.putAtt(ncid,id_ff,'short_name','wspeed');
netcdf.putAtt(ncid,id_ff,'units','m s-1');
netcdf.putAtt(ncid,id_ff,'valid_range',[nanmin(data.ff(:)), nanmax(data.ff(:))]);
netcdf.putAtt(ncid,id_ff,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_ff,'WIND.SPEED.SURFACE_SOURCE','Vaisala weather station WXT520 or WXT530');
netcdf.putAtt(ncid,id_ff,'comment',['Speed is the magnitude of velocity. Wind '...
                                    'is defined as a two-dimensional (horizontal) '...
                                    'air velocity vector, with no vertical '...
                                    'component. (Vertical motion in the atmosphere '...
                                    'has the standard name upward_air_velocity.) '...
                                    'The wind speed is the magnitude of the wind velocity.']);

id_fff = netcdf.defVar(ncid,'wdir','nc_float',did_time);
netcdf.putAtt(ncid,id_fff,'GEOMS_name','WIND.DIRECTION.SURFACE');
netcdf.putAtt(ncid,id_fff,'long_name','wind_from_direction');
netcdf.putAtt(ncid,id_fff,'rpg_name','fff');
netcdf.putAtt(ncid,id_fff,'standard_name','WIND.DIRECTION measured by the meteo-station attached to the radar');
netcdf.putAtt(ncid,id_fff,'short_name','wdir');
netcdf.putAtt(ncid,id_fff,'units','degrees');
netcdf.putAtt(ncid,id_fff,'valid_range',[nanmin(data.fff(:)), nanmax(data.fff(:))]);
netcdf.putAtt(ncid,id_fff,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_fff,'wind.direction.surface_source',['Vaisala weather '...
                           'station WXT520 or WXT530']);
netcdf.putAtt(ncid,id_fff,'comment',['Wind is defined as a two-dimensional '...
                                     '(horizontal) air velocity vector, with '...
                                     'no vertical component. (Vertical motion '...
                                     'in the atmosphere has the standard name '...
                                     'upward_air_velocity.) In meteorological '...
                                     'reports, the direction of the wind vector '...
                                     'is usually (but not always) given as the '...
                                     'direction from which it is blowing '...
                                     '(wind_from_direction) (westerly, northerly, etc.). '...
                                     'In other contexts, such as atmospheric '...
                                     'modelling, it is often natural to give '...
                                     'the direction in the usual manner of '...
                                     'vectors as the heading or the direction '...
                                     'to which it is blowing (wind_to_direction) '...
                                     '(eastward, southward, etc.) "from_direction" '...
                                     'is used in the construction X_from_direction '...
                                     'and indicates the direction from which '...
                                     'the velocity vector of X is coming.']);
                       
id_Tb = netcdf.defVar(ncid,'tb','nc_float',did_time);
netcdf.putAtt(ncid,id_Tb,'GEOMS_name','TEMPERATURE.BRIGHTNESS');
netcdf.putAtt(ncid,id_Tb,'long_name','brightness_temperature');
netcdf.putAtt(ncid,id_Tb,'rpg_name','Tb');
netcdf.putAtt(ncid,id_Tb,'standard_name',['Brightness temperature at 89-GHz detected '...
                                      'by the passive chanal of the receiver '...
                                      'antenna of the radar']);
netcdf.putAtt(ncid,id_Tb,'short_name','tb');
netcdf.putAtt(ncid,id_Tb,'units','K');
netcdf.putAtt(ncid,id_Tb,'valid_range',[nanmin(data.Tb(:)), nanmax(data.Tb(:))]);
netcdf.putAtt(ncid,id_Tb,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_Tb,'comment',['The brightness temperature of a body is '...
                                    'the temperature of a black body which '...
                                    'radiates the same power per unit solid '...
                                    'angle per unit area. Brightness Temperature '...
                                    'measurements at 89-GHz by the radar']);

id_lwp = netcdf.defVar(ncid,'lwp','nc_float',did_time);
netcdf.putAtt(ncid,id_lwp,'GEOMS_name','LIQUID.WATER.PATH');
netcdf.putAtt(ncid,id_lwp,'long_name','liquid_water_path');
netcdf.putAtt(ncid,id_lwp,'rpg_name','lwp');
netcdf.putAtt(ncid,id_lwp,'standard_name','Liquid water path (lwp) calculated by RPG software');
netcdf.putAtt(ncid,id_lwp,'short_name','lwp');
netcdf.putAtt(ncid,id_lwp,'units','g m-2');
netcdf.putAtt(ncid,id_lwp,'valid_range',[nanmin(data.lwp(:)), nanmax(data.lwp(:))]);
netcdf.putAtt(ncid,id_lwp,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_lwp,'comment',['A measure of the weight of the liquid '...
                                     'water droplets in the atmosphere above '...
                                     'a unit surface area on the earth. The '...
                                     'liquid  water path is calculated from '...
                                     'the tb measurement of the 89-GHz chanal. '...
                                     'The retrieval was developed by RPG and '...
                                     'is based on a nural network approach']);

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
netcdf.putAtt(ncid,id_T_trans,'GEOMS_name','TRANSMITTER.TEMPERATUR');
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
netcdf.putAtt(ncid,id_T_rec,'GEOMS_name','RECEIVER.TEMPERATUR');
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
netcdf.putAtt(ncid,id_T_pc,'GEOMS_name','COMPUTER.TEMPERATUR');
netcdf.putAtt(ncid,id_T_pc,'long_name','platform_computer_temperature');
netcdf.putAtt(ncid,id_T_pc,'rpg_name','T_pc');
netcdf.putAtt(ncid,id_T_pc,'standard_name','Temperature of the radar controlling PC');
netcdf.putAtt(ncid,id_T_pc,'short_name','t_pc');
netcdf.putAtt(ncid,id_T_pc,'units','K');
netcdf.putAtt(ncid,id_T_pc,'valid_range',[nanmin(data.T_pc(:)), nanmax(data.T_pc(:))]);
netcdf.putAtt(ncid,id_T_pc,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_T_pc,'comment',['The radar controllin PC temperature '...
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
                                 
%%%%%%%% multi-D variables

id_Ze = netcdf.defVar(ncid,'ze','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_Ze,'GEOMS_name','RADAR.REFLECTIVITY.FACTOR_VV');
netcdf.putAtt(ncid,id_Ze,'long_name','equivalent_reflectivity_factor');
netcdf.putAtt(ncid,id_Ze,'rpg_name','Ze');
netcdf.putAtt(ncid,id_Ze,'standard_name','Equivalent radar reflectivity factor at vertical polarisation');
netcdf.putAtt(ncid,id_Ze,'units','mm6 m-3');
netcdf.putAtt(ncid,id_Ze,'valid_range',[nanmin(data.Ze(:)), nanmax(data.Ze(:))]);
netcdf.putAtt(ncid,id_Ze,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_Ze,'comment',['"Equivalent reflectivity factor" is the '...
                                    'radar reflectivity factor that is calculated '...
                                    'from the measured radar return power '...
                                    'assuming the target is composed of liquid '...
                                    'water droplets whose diameter is less than '...
                                    'one tenth of the radar wavelength, i.e., '...
                                    'treating the droplets as Rayleigh scatterers. '...
                                    'The actual radar reflectivity factor would '...
                                    'depend on the size distribution and '...
                                    'composition of the particles within the '...
                                    'target volume and these are often unknown.']);
if isfield(data, 'Ze_label') % Ze corrected, adding note
    netcdf.putAtt(ncid,id_Ze,'comment2',data.Ze_label);
    netcdf.putAtt(ncid,id_Ze,'corretion_dB',data.Ze_corr);
end

id_vm = netcdf.defVar(ncid,'vm','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_vm,'GEOMS_name','DOPPLER.VELOCITY_MEAN');
netcdf.putAtt(ncid,id_vm,'long_name','mean_doppler_velocity');
netcdf.putAtt(ncid,id_vm,'rpg_name','vm');
netcdf.putAtt(ncid,id_vm,'standard_name','Mean Doppler velocity at vertical polarisation');
netcdf.putAtt(ncid,id_vm,'short_name','vm');
netcdf.putAtt(ncid,id_vm,'units','m s-1');
netcdf.putAtt(ncid,id_vm,'valid_range',[nanmin(data.vm(:)), nanmax(data.vm(:))]);
netcdf.putAtt(ncid,id_vm,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_vm,'comment',['The Doppler velocity is the radial component of '...
                                    'the velocity vector of a scattering object as '...
                                    'observed by a remote sensor, such as a Doppler. '...
                                    'NOTE: Negative velocities indicate particles '...
                                    'motion towards the radar'])

id_sigma = netcdf.defVar(ncid,'sw','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_sigma,'GEOMS_name','DOPPLER.SPECTRUM_WIDTH');
netcdf.putAtt(ncid,id_sigma,'long_name','doppler_spectrum_width');
netcdf.putAtt(ncid,id_sigma,'rpg_name','sigma');
netcdf.putAtt(ncid,id_sigma,'standard_name','Doppler spectrum width at vertical polarisation');
netcdf.putAtt(ncid,id_sigma,'short_name','sw');
netcdf.putAtt(ncid,id_sigma,'units','m s-1');
netcdf.putAtt(ncid,id_sigma,'valid_range',[nanmin(data.sigma(:)), nanmax(data.sigma(:))]);
netcdf.putAtt(ncid,id_sigma,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_sigma,'comment',['Doppler scectrum width of the measured Doppler spectrum'...
                                       'at vertical polarisation.'])

id_skew = netcdf.defVar(ncid,'skew','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_skew,'GEOMS_name','DOPPLER.SPECTRUM_SKEWNESS');
netcdf.putAtt(ncid,id_skew,'long_name','doppler_spectrum_skewnes');
netcdf.putAtt(ncid,id_skew,'rpg_name','skew');
netcdf.putAtt(ncid,id_skew,'standard_name','Doppler spectrum skewness at vertical polarisatio');
netcdf.putAtt(ncid,id_skew,'short_name','skew');
netcdf.putAtt(ncid,id_skew,'valid_range',[nanmin(data.skew(:)), nanmax(data.skew(:))]);
netcdf.putAtt(ncid,id_skew,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_skew,'comment',['Doppler spectrum skewness of measued Doppler spectrum '...
                                      'at vertical polarisation.'])
%Included by Bravo-Aranda, J.A. JABA
if data.DualPol > 0
    id_Ze_hv = netcdf.defVar(ncid,'ze_hv','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_Ze_hv,'GEOMS_name','RADAR.REFLECTIVITY.FACTOR_HV');
    netcdf.putAtt(ncid,id_Ze_hv,'long_name','equivalent_reflectivity_factor');
    netcdf.putAtt(ncid,id_Ze_hv,'rpg_name','Ze_hv');
    netcdf.putAtt(ncid,id_Ze_hv,'standard_name','Equivalent radar reflectivity factor of the cross polarized chanal');
    netcdf.putAtt(ncid,id_Ze_hv,'units','mm6 m-3');
    netcdf.putAtt(ncid,id_Ze_hv,'valid_range',[nanmin(data.Ze(:)), nanmax(data.Ze(:))]);
    netcdf.putAtt(ncid,id_Ze_hv,'fill_value',str_fill_value);
    netcdf.putAtt(ncid,id_Ze_hv,'comment',['"Equivalent reflectivity factor" is the '...
                                        'radar reflectivity factor that is calculated '...
                                        'from the measured radar return power '...
                                        'assuming the target is composed of liquid '...
                                        'water droplets whose diameter is less than '...
                                        'one tenth of the radar wavelength, i.e., '...
                                        'treating the droplets as Rayleigh scatterers. '...
                                        'The actual radar reflectivity factor would '...
                                        'depend on the size distribution and '...
                                        'composition of the particles within the '...
                                        'target volume and these are often unknown.']);
    
    id_ldr = netcdf.defVar(ncid,'ldr','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_ldr,'GEOMS_name','LINEAR.DEPOLARIZATION.RATIO');
    netcdf.putAtt(ncid,id_ldr,'long_name','linear_depolarisation_ratio');
    netcdf.putAtt(ncid,id_ldr,'rpg_name','ldr');
    netcdf.putAtt(ncid,id_ldr,'standard_name','Linear depolarization ratio');
    netcdf.putAtt(ncid,id_ldr,'short_name','ldr');
    netcdf.putAtt(ncid,id_ldr,'unite','mm6 m-3');
    netcdf.putAtt(ncid,id_ldr,'valid_range',[nanmin(data.LDR(:)), nanmax(data.LDR(:))]);
    netcdf.putAtt(ncid,id_ldr,'fill values',str_fill_value);
    netcdf.putAtt(ncid,id_ldr,'comment',['The ratio of the power received in '...
                                         'the orthogonal, or cross-polarized, '...
                                         'channel to that received in the '...
                                         'transmission, or copolarized, channel'...
                                         'of a dual- channel radar, when a '...
                                         'linearly polarized signal is transmitted. '...
                                         'L_dr = Ze_hv/Ze_vv']);
end
if data.DualPol == 2
    id_xcorr = netcdf.defVar(ncid,'rho_hv','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_xcorr,'GEOMS_name','CORRELATION.COEFFICIENT.');
    netcdf.putAtt(ncid,id_xcorr,'long_name','correlation_coefficient');
    netcdf.putAtt(ncid,id_xcorr,'rpg_name','CorrCoeff');    
    netcdf.putAtt(ncid,id_xcorr,'standard_name','co-cross-channel correlation coefficient');
    netcdf.putAtt(ncid,id_xcorr,'short_name','rho_hv');        
    %netcdf.putAtt(ncid,id_xcorr,'valid_range',[min(data.xcorr(:)), max(data.xcorr(:))]);
    netcdf.putAtt(ncid,id_xcorr,'fill_value',str_fill_value);
    
    id_difphase = netcdf.defVar(ncid,'phi_dp','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_difphase,'GEOMS_name','DIFFERENTIAL.PHASE');
    netcdf.putAtt(ncid,id_difphase,'long_name','differential_phase');
    netcdf.putAtt(ncid,id_difphase,'rpg_name','difphase');    
    netcdf.putAtt(ncid,id_difphase,'standard_name','co-cross-channel differential phase');   
    netcdf.putAtt(ncid,id_difphase,'short_name','phi_dp');        
    netcdf.putAtt(ncid,id_difphase,'unite','degree');   
    %netcdf.putAtt(ncid,id_difphase,'valid_range',[min(data.difphase(:)), max(data.difphase(:))]); 
    netcdf.putAtt(ncid,id_difphase,'fill_value',str_fill_value); 
end

id_QualFlag = netcdf.defVar(ncid,'pro_qf','nc_byte',[did_range,did_time]);
netcdf.putAtt(ncid,id_QualFlag,'GEOMS_name','FLAG.PROCESSING.QUALITY');
netcdf.putAtt(ncid,id_QualFlag,'long_name','quality_falg_data_processing');
netcdf.putAtt(ncid,id_QualFlag,'standard_name','Quality flag, added in the additional data processing to alert for known issues');
netcdf.putAtt(ncid,id_QualFlag,'short_name','pro_qf');
netcdf.putAtt(ncid,id_QualFlag,'comment', ...
    ['This variable contains information on anything that might impact the quality ', ...
    'of the data at each pixel. Must be converted into three bit binary string. ', ...
	'If 0, i.e. dec2bin(QualityFlag,3) = 000, none of the included issues were ', ...
    'found. The definitions of each bit are given in the definition attribute.']);
netcdf.putAtt(ncid,id_QualFlag,'definition', ...
    ['If 2^0 bit is 1: this range gate is known to have aritifical spikes occurring ', ...
     'If 2^1 bit is 1: aircraft or other known flying non-meteorological object ', ...
     'If 2^2 bit is 1: wet-radome']);

id_Aliasmask = netcdf.defVar(ncid,'alias_mask','nc_byte',[did_range,did_time]);
netcdf.putAtt(ncid,id_Aliasmask,'GEOMS_name','MASK.PROCESSING.ALIAISING');
netcdf.putAtt(ncid,id_Aliasmask,'long_name','alaising_mask');
netcdf.putAtt(ncid,id_Aliasmask,'standard_name','Mask array indicating in which bin aliasing occurs');
netcdf.putAtt(ncid,id_Aliasmask,'comment',['AliasMask indicates the bins where aliasing is detected. '...
                                           '0 = no aliasing; 1 = aliasing occurs. Only applicable if '...
                                           'variable AntiAlias is 2.']);

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
netcdf.defVarDeflate(ncid,id_incel,true,true,9);
netcdf.defVarDeflate(ncid,id_incea,true,true,9);

%radar moments
netcdf.defVarDeflate(ncid,id_Ze,true,true,9);
netcdf.defVarDeflate(ncid,id_vm,true,true,9);
netcdf.defVarDeflate(ncid,id_sigma,true,true,9);
netcdf.defVarDeflate(ncid,id_skew,true,true,9);
if data.DualPol > 0
    netcdf.defVarDeflate(ncid,id_ldr,true,true,9); %JABA    
    netcdf.defVarDeflate(ncid,id_Ze_hv,true,true,9); %LP
end 
if data.DualPol == 2
    netcdf.defVarDeflate(ncid,id_difphase,true,true,9); %JABA    
    netcdf.defVarDeflate(ncid,id_xcorr,true,true,9); %JABA  
end
%processing related variables 
netcdf.defVarDeflate(ncid,id_QualFlag,true,true,9);
netcdf.defVarDeflate(ncid,id_Aliasmask,true,true,9);

netcdf.endDef(ncid);

%% ####################### put variables into file

% yourmatrix(isnan(yourmatrix)) = somevalue;

% scalars
netcdf.putVar(ncid,id_lat,0,data.Lat);
netcdf.putVar(ncid,id_lon,0,data.Lon);
netcdf.putVar(ncid,id_MSL,0,data.MSL);
netcdf.putVar(ncid,id_freq,0, 299792458/(data.freq * 1e9) ); 
netcdf.putVar(ncid,id_wl,0,data.freq * 1e9);
netcdf.putVar(ncid,id_HPBW,0,data.HPBW);

% range dependet
netcdf.putVar(ncid,id_range,0,data.n_levels,data.range);

% chrip seq dependent variables
netcdf.putVar(ncid,id_range_offsets,0,data.no_chirp_seq,data.range_offsets);
netcdf.putVar(ncid,id_DoppMax,0,data.no_chirp_seq,data.DoppMax);
% netcdf.putVar(ncid,id_AntiAlias,0,data.AntiAlias);

% time dependent variables
% netcdf.putVar(ncid,id_sampleTms,0,data.totsamp,data.sampleTms);
netcdf.putVar(ncid,id_time,0,data.totsamp,data.time);
temporary_data = data.RR ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_RR,0,data.totsamp,temporary_data);
clear temporary_data ;
temporary_data = data.rh ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_rh,0,data.totsamp,temporary_data);
clear temporary_data ;
temporary_data = data.T_env ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_T_env,0,data.totsamp, temporary_data);
clear temporary_data ;
temporary_data = data.pres ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_pres,0,data.totsamp, temporary_data);
clear temporary_data ;
temporary_data = data.ff * (1000/3600) ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_ff,0,data.totsamp,temporary_data);
clear temporary_data ;
temporary_data = data.fff ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_fff,0,data.totsamp,temporary_data);
clear temporary_data ;
temporary_data = data.Tb ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_Tb,0,data.totsamp,temporary_data);
clear temporary_data ;
temporary_data = data.lwp ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_lwp,0,data.totsamp,temporary_data);
netcdf.putVar(ncid,id_status,0,data.totsamp,data.status);
clear temporary_data ;
temporary_data = data.TransPow ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_TransPow,0,data.totsamp,temporary_data);
clear temporary_data ;
temporary_data = data.T_trans ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_T_trans,0,data.totsamp,temporary_data);
clear temporary_data ;
temporary_data = data.T_rec ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_T_rec,0,data.totsamp,temporary_data);
clear temporary_data ;
temporary_data = data.T_pc ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_T_pc,0,data.totsamp,temporary_data);
netcdf.putVar(ncid,id_QF,0,data.totsamp,data.QF);

% multidimensional variables
clear temporary_data ;
temporary_data = data.Ze ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_Ze,[0,0],[data.n_levels,data.totsamp],temporary_data');
clear temporary_data ;
temporary_data = data.vm ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_vm,[0,0],[data.n_levels,data.totsamp],temporary_data');
clear temporary_data ;
temporary_data = data.sigma ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_sigma,[0,0],[data.n_levels,data.totsamp],temporary_data');
clear temporary_data ;
temporary_data = data.skew ;
temporary_data(isnan(temporary_data))=fill_value ;
netcdf.putVar(ncid,id_skew,[0,0],[data.n_levels,data.totsamp],temporary_data');
netcdf.putVar(ncid,id_QualFlag,[0,0],[data.n_levels,data.totsamp],data.QualFlag');
netcdf.putVar(ncid,id_Aliasmask,[0,0],[data.n_levels,data.totsamp],data.Aliasmask');
netcdf.putVar(ncid,id_incel,0,data.totsamp,data.reserved(:,2));
netcdf.putVar(ncid,id_incea,0,data.totsamp,data.reserved(:,3));

if data.DualPol > 0
    clear temporary_data ;
    temporary_data = data.LDR ;
    temporary_data(isnan(temporary_data))=fill_value ;
    netcdf.putVar(ncid,id_ldr,[0,0],[data.n_levels,data.totsamp],temporary_data'); %JABA
    clear temporary_data ;
    temporary_data = data.Ze_hv ;
    temporary_data(isnan(temporary_data))=fill_value ;
    netcdf.putVar(ncid,id_Ze_hv,[0,0],[data.n_levels,data.totsamp],temporary_data'); %LP
end

if data.DualPol == 2
    disp('WARNING! skipping writing id_xcorr and difphase into files, variables not available (yet??)')
    % netcdf.putVar(ncid,id_xcorr,[0,0],[data.n_levels,data.totsamp],data.xcorr'); %JABA
    % netcdf.putVar(ncid,id_difphase,[0,0],[data.n_levels,data.totsamp],data.difphase'); %JABA
end

netcdf.close(ncid);

end
