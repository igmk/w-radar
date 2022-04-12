function write_joyrad94_data_2_nc_geoms(data, outfile, config)

% this function writes joyrad94 data into netcdf4
% Changes of all the 

%% ################## Create a netCDF file.

time_start = datestr(datenum([2001 1 1 0 0 0]) + data.time(1)/(86400),'yyyymmddTHHMMSS') ;
time_end   = datestr(datenum([2001 1 1 0 0 0]) + data.time(end)/(86400),'yyyymmddTHHMMSS') ;

%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified by A.Pschera 2021-11-08 - added especially '.nc'
outfile_geoms = ([config.outputpath_tree,'/','groundbased_','radar_profiler_',config.nickradar,'_',config.file_association,'_',config.nickstation,'_',time_start,'_',time_end,'_',config.data_version_no,'.nc']);
%%%%%%%%%%%%%%%%%%%%%%%%%

ncid = netcdf.create(outfile_geoms,'NETCDF4');


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
netcdf.putAtt(ncid,glob,'DATA_PROCESSOR',config.data_processor);
netcdf.putAtt(ncid,glob,'DATA_START_DATE',time_start);
netcdf.putAtt(ncid,glob,'DATA_STOP_DATE',time_end);
netcdf.putAtt(ncid,glob,'DATA_FILE_VERSION',config.data_version_no);
netcdf.putAtt(ncid,glob,'DATA_MODIFICATIONS',config.data_modification);
netcdf.putAtt(ncid,glob,'DATA_CAVEATS',config.data_caveats);
netcdf.putAtt(ncid,glob,'DATA_RULES_OF_USE',config.data_rules_of_use);
netcdf.putAtt(ncid,glob,'DATA_ACKNOWLEDGEMENT',config.data_acknowledgement);
netcdf.putAtt(ncid,glob,'DATA_QUALITY',config.data_quality);
netcdf.putAtt(ncid,glob,'DATA_TEMPLATE',config.data_tamplate);

netcdf.putAtt(ncid,glob,'FILE_NAME',['groundbased_',...
                                     'radar_profiler','_',...
                                     config.nickradar,'_',...
                                     config.file_association,'_',...
                                     config.nickstation,'_',...
                                     time_start,'_',time_end,'_',...
                                     config.data_version_no]);
netcdf.putAtt(ncid,glob,'FILE_GENERATION_DATE',datestr(now,'yyyymmddTHHMMSS'));
netcdf.putAtt(ncid,glob,'FILE_ACCESS',config.file_access);
netcdf.putAtt(ncid,glob,'FILE_PROJECT_ID',config.file_project_id);
netcdf.putAtt(ncid,glob,'FILE_DOI',config.file_doi);
netcdf.putAtt(ncid,glob,'FILE_ASSOCIATION',config.file_association);
netcdf.putAtt(ncid,glob,'FILE_META_VERSION',config.file_metadata_version);

netcdf.putAtt(ncid,glob,'FILL_VALUE',str_fill_value);
netcdf.putAtt(ncid,glob,'INSTRUMENT_MODEL',model);
netcdf.putAtt(ncid,glob,'MDF_PROGRAM_USED',data.progname);


%% ################ get variable ids and add attributes

%%%%%%%%%% scalar variables

id_lat = netcdf.defVar(ncid,'LATITUDE.INSTRUMENT','nc_float',did_scalar);
netcdf.putAtt(ncid,id_lat,'VAR_NAME','LATITUDE.INSTRUMENT');
netcdf.putAtt(ncid,id_lat,'VAR_DESCRIPTION',['Inst. geolocation. Latitude '...
                                             'north (decimal degrees) of the '...
                                             'location of the instrument ' ...
                                             '(+ for north; - for south)']);
netcdf.putAtt(ncid,id_lat,'VAR_NOTES',' ');
netcdf.putAtt(ncid,id_lat,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_lat,'VAR_DEPEND','CONSTANT');
netcdf.putAtt(ncid,id_lat,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_lat,'VAR_UNITS','deg');
netcdf.putAtt(ncid,id_lat,'VAR_SI_CONVERSION','0.0;1.74533E-2;rad');
netcdf.putAtt(ncid,id_lat,'VAR_VALIS_MIN','-90.0');
netcdf.putAtt(ncid,id_lat,'VAR_VALIS_MAX','90.0');
netcdf.putAtt(ncid,id_lat,'VAR_FILL_VALUE','LATITUDE.INSTRUMENT');
netcdf.putAtt(ncid,id_lat,'long_name','latitude');
netcdf.putAtt(ncid,id_lat,'short_name','lat');

id_lon = netcdf.defVar(ncid,'LONGITUDE.INSTRUMENT','nc_float',did_scalar);
netcdf.putAtt(ncid,id_lon,'VAR_NAME','LONGITUDE.INSTRUMENT');
netcdf.putAtt(ncid,id_lon,'VAR_DESCRIPTION',['Inst. geolocation. Longitude '...
                                             'east (decimal degrees) of the ' ...
                                             'location of the instrument '...
                                             '(+ for east; - for west)']);
netcdf.putAtt(ncid,id_lon,'VAR_NOTES',' ');
netcdf.putAtt(ncid,id_lon,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_lon,'VAR_DEPEND','CONSTANT');
netcdf.putAtt(ncid,id_lon,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_lon,'VAR_UNITS','deg');
netcdf.putAtt(ncid,id_lon,'VAR_SI_CONVERSION','0.0;1.74533E-2;rad');
netcdf.putAtt(ncid,id_lon,'VAR_VALIS_MIN','-180.0');
netcdf.putAtt(ncid,id_lon,'VAR_VALIS_MAX','180.0');
netcdf.putAtt(ncid,id_lon,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_lon,'long_name','longitude');
netcdf.putAtt(ncid,id_lon,'short_name','lon');

id_MSL = netcdf.defVar(ncid,'ALTITUDE.INSTRUMENT','nc_float',did_scalar);
netcdf.putAtt(ncid,id_MSL,'VAR_NAME','ALTITUDE.INSTRUMENT');
netcdf.putAtt(ncid,id_MSL,'VAR_DESCRIPTION',['Inst. geolocation. Altitude of '...
                                             'the instrument relative to the '...
                                             'location site']);
netcdf.putAtt(ncid,id_MSL,'VAR_NOTES',' ');
netcdf.putAtt(ncid,id_MSL,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_MSL,'VAR_DEPEND','CONSTANT');
netcdf.putAtt(ncid,id_MSL,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_MSL,'VAR_UNITS','m');
netcdf.putAtt(ncid,id_MSL,'VAR_SI_CONVERSION','0.0;1.0;m');
netcdf.putAtt(ncid,id_MSL,'VAR_VALIS_MIN','0.0');
netcdf.putAtt(ncid,id_MSL,'VAR_VALIS_MAX','10000.0');
netcdf.putAtt(ncid,id_MSL,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_MSL,'long_name','altitude');
netcdf.putAtt(ncid,id_MSL,'short_name','zsl');

id_freq = netcdf.defVar(ncid,'FREQUENCY','nc_float',did_scalar);
netcdf.putAtt(ncid,id_freq,'VAR_NAME','FREQUENCY');
netcdf.putAtt(ncid,id_freq,'VAR_DESCRIPTION',['Central transmission frequency '...
                                              'of the electromagnetic wave '...
                                              'emitted by radar']);
netcdf.putAtt(ncid,id_freq,'VAR_NOTES',['Frequency can be converted into ' ...
                                        'WAVELENGTH = c/FREQUENCY; where ' ...
                                        'c is the speed of light']);
netcdf.putAtt(ncid,id_freq,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_freq,'VAR_DEPEND','CONSTANT');
netcdf.putAtt(ncid,id_freq,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_freq,'VAR_UNITS','m-1');
netcdf.putAtt(ncid,id_freq,'VAR_SI_CONVERSION','0.0;1.0;m');
netcdf.putAtt(ncid,id_freq,'VAR_VALIS_MIN','0.0');
netcdf.putAtt(ncid,id_freq,'VAR_VALIS_MAX','0.01');
netcdf.putAtt(ncid,id_freq,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_freq,'long_name','radiation_frequency');
netcdf.putAtt(ncid,id_freq,'short_name','freq');

id_wl = netcdf.defVar(ncid,'WAVELENGTH','nc_float',did_scalar);
netcdf.putAtt(ncid,id_wl,'VAR_NAME','WAVELENGTH');
netcdf.putAtt(ncid,id_wl,'VAR_DESCRIPTION',['Information about operating radar '...
                                           'wavelengh was calculated based on '...
                                           'the frequency information']);
netcdf.putAtt(ncid,id_wl,'VAR_NOTES',['Frequency can be converted into ' ...
                                        'WAVELENGTH = c/FREQUENCY; where ' ...
                                        'c is the speed of light']);
netcdf.putAtt(ncid,id_wl,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_wl,'VAR_DEPEND','CONSTANT');
netcdf.putAtt(ncid,id_wl,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_wl,'VAR_UNITS','m');
netcdf.putAtt(ncid,id_wl,'VAR_SI_CONVERSION','0.0;1.0;m');
netcdf.putAtt(ncid,id_wl,'VAR_VALIS_MIN','9.0+e+10f');
netcdf.putAtt(ncid,id_wl,'VAR_VALIS_MAX','9.5e+10f');
netcdf.putAtt(ncid,id_wl,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_wl,'long_name','radiation_wavelength');
netcdf.putAtt(ncid,id_wl,'short_name','wl');


id_HPBW = netcdf.defVar(ncid,'BEAM.WIDTH','nc_float',did_scalar);
netcdf.putAtt(ncid,id_HPBW,'VAR_NAME','BEAM.WIDTH');
netcdf.putAtt(ncid,id_HPBW,'VAR_DESCRIPTION',['Antenna half power beam width '...
                                              '(re. https://www.radartutorial.eu'...
                                              '/06.antennas/an08.en.html)']);
netcdf.putAtt(ncid,id_HPBW,'VAR_NOTES',['Half power beam width is the angle between '...
                                      'the half-power (-3 dB) points of the main '...
                                      'lobe, when referenced to the peak effective '...
                                      'radiated power of the main lobe']);
netcdf.putAtt(ncid,id_HPBW,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_HPBW,'VAR_DEPEND','CONSTANT');
netcdf.putAtt(ncid,id_HPBW,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_HPBW,'VAR_UNITS','deg');
netcdf.putAtt(ncid,id_HPBW,'VAR_SI_CONVERSION','0.0;1.74533E-2;rad');
netcdf.putAtt(ncid,id_HPBW,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_HPBW,'VAR_VALIS_MAX','2.0f');
netcdf.putAtt(ncid,id_HPBW,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_HPBW,'long_name','antenna_beam_width');
netcdf.putAtt(ncid,id_HPBW,'short_name','hpbw');

id_AntiAlias = netcdf.defVar(ncid,'ALIAS_FLAG','nc_byte',did_scalar);
netcdf.putAtt(ncid,id_AntiAlias,'VAR_NAME','ALIAS_FLAG');
netcdf.putAtt(ncid,id_AntiAlias,'VAR_DESCRIPTION',['The flag index shows: '...
                                                   '0 = no dealiasing applied, '... 
                                                   '1 = dealiasing by RPG, '...
                                                   '2 = dealiasing by the applied '...
                                                   'code (see DATA_SOURCE)']);
netcdf.putAtt(ncid,id_AntiAlias,'VAR_NOTES','');
netcdf.putAtt(ncid,id_AntiAlias,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_AntiAlias,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_AntiAlias,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_AntiAlias,'VAR_UNITS','1');
netcdf.putAtt(ncid,id_AntiAlias,'VAR_SI_CONVERSION','0.0;0.01;');
netcdf.putAtt(ncid,id_AntiAlias,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_AntiAlias,'VAR_VALIS_MAX','2.0f');
netcdf.putAtt(ncid,id_AntiAlias,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_AntiAlias,'long_name','Quality flag for dealiasing');

                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% range variables %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len_range = num2str(length(data.range)) ;
no_chirps = num2str(length(data.no_chirp_seq)) ;

id_range = netcdf.defVar(ncid,'RANGE','nc_float',did_range);
netcdf.putAtt(ncid,id_range,'VAR_NAME','RANGE');
netcdf.putAtt(ncid,id_range,'VAR_DESCRIPTION',['Distance from the Sensor to '...
                                               'Center of each Range Gate '...
                                               'along the line of sight']);
netcdf.putAtt(ncid,id_range,'VAR_NOTES','');
netcdf.putAtt(ncid,id_range,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_range,'VAR_DEPEND','RANGE');
netcdf.putAtt(ncid,id_range,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_range,'VAR_UNITS','m');
netcdf.putAtt(ncid,id_range,'VAR_SI_CONVERSION','0.0;1.0;m');
netcdf.putAtt(ncid,id_range,'VAR_VALIS_MIN','0.0');
netcdf.putAtt(ncid,id_range,'VAR_VALIS_MAX','20000.0');
netcdf.putAtt(ncid,id_range,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_range,'long_name','range');
netcdf.putAtt(ncid,id_range,'short_name','range');

%%%%%%%% chirp_seq_dependent variables

id_DoppMax = netcdf.defVar(ncid,'NYQUIST.VELOCITY','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_DoppMax,'VAR_NAME','NYQUIST.VELOCITY');
netcdf.putAtt(ncid,id_DoppMax,'VAR_DESCRIPTION',['The Nyquist velocity is the '...
                                               'maximum velocity that can be '...
                                               'correctly displayed by a '...
                                               'Doppler radar, and this is '...
                                               'dependent on the wavelength '...
                                               'and frequency of eleoctromagnetic '...
                                               'wave emitted by the radar']);
netcdf.putAtt(ncid,id_DoppMax,'VAR_NOTES','');
netcdf.putAtt(ncid,id_DoppMax,'VAR_SIZE',no_chirps);
netcdf.putAtt(ncid,id_DoppMax,'VAR_DEPEND','CHIRP.SEQUENCES.INDEX');
netcdf.putAtt(ncid,id_DoppMax,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_DoppMax,'VAR_UNITS','m s-1');
netcdf.putAtt(ncid,id_DoppMax,'VAR_SI_CONVERSION','0.0;1.0;m s-1');
netcdf.putAtt(ncid,id_DoppMax,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_DoppMax,'VAR_VALIS_MAX','50.0f');
netcdf.putAtt(ncid,id_DoppMax,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_DoppMax,'long_name','nyquist_velocity');
netcdf.putAtt(ncid,id_DoppMax,'short_name','nqv');


id_range_offsets = netcdf.defVar(ncid,'RANGE_OFFSET','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_range_offsets,'VAR_NAME','RANGE_OFFSET');
netcdf.putAtt(ncid,id_range_offsets,'VAR_DESCRIPTION',['Range offsets between '...
                                                       'different chirps when '...
                                                       'measuring with FMCW-radars.']);
netcdf.putAtt(ncid,id_range_offsets,'VAR_NOTES',['The range(range_offsets) '...
                                                'will give you the range where a '...
                                                'new chirp sequence starts. '...
                                                'range_offsets counts from '...
                                                '1 to n_levels.']);
netcdf.putAtt(ncid,id_range_offsets,'VAR_SIZE',no_chirps);
netcdf.putAtt(ncid,id_range_offsets,'VAR_DEPEND','CHIRP.SEQUENCES.INDEX');
netcdf.putAtt(ncid,id_range_offsets,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_range_offsets,'VAR_UNITS','m');
netcdf.putAtt(ncid,id_range_offsets,'VAR_SI_CONVERSION','0.0;1.0;m s-1');
netcdf.putAtt(ncid,id_range_offsets,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_range_offsets,'VAR_VALIS_MAX','100.0f');
netcdf.putAtt(ncid,id_range_offsets,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_range_offsets,'long_name','range_offset');
netcdf.putAtt(ncid,id_range_offsets,'short_name','range_offset');
                     

id_vm_res = netcdf.defVar(ncid,'DOPPLER.VELOCITY_RESOLUTION','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_vm_res,'VAR_NAME','DOPPLER.VELOCITY_RESOLUTION');
netcdf.putAtt(ncid,id_vm_res,'VAR_DESCRIPTION',['The resolution of the Doppler '...
                                                'spectrum between +/- Nyquist '...
                                                'velocity for each chirp']);
netcdf.putAtt(ncid,id_vm_res,'VAR_NOTES','');
netcdf.putAtt(ncid,id_vm_res,'VAR_SIZE',no_chirps);
netcdf.putAtt(ncid,id_vm_res,'VAR_DEPEND','CHIRP.SEQUENCES.INDEX');
netcdf.putAtt(ncid,id_vm_res,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_vm_res,'VAR_UNITS','m s-1');
netcdf.putAtt(ncid,id_vm_res,'VAR_SI_CONVERSION','0.0;1.0;m s-1');
netcdf.putAtt(ncid,id_vm_res,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_vm_res,'VAR_VALIS_MAX','50.0f');
netcdf.putAtt(ncid,id_vm_res,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_vm_res,'long_name','doppler_velocity_resolution');
netcdf.putAtt(ncid,id_vm_res,'short_name','svm_res');

id_int_time = netcdf.defVar(ncid,'INTEGRATION.TIME','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_int_time,'VAR_NAME','INTEGRATION.TIME');
netcdf.putAtt(ncid,id_int_time,'VAR_DESCRIPTION',['Duration or integration ' ...
                                                  'time per single chirp.']);
netcdf.putAtt(ncid,id_int_time,'VAR_NOTES',['The duration or integration time '...
                                            'might differ from the single '...
                                            'time steps in DATETIME. This is '...
                                            'some hard or software also needs ' ...
                                            'time to finish processes.']);
netcdf.putAtt(ncid,id_int_time,'VAR_SIZE',no_chirps);
netcdf.putAtt(ncid,id_int_time,'VAR_DEPEND','CHIRP.SEQUENCES.INDEX');
netcdf.putAtt(ncid,id_int_time,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_int_time,'VAR_UNITS','s');
netcdf.putAtt(ncid,id_int_time,'VAR_SI_CONVERSION','0.0;86400.0;s');
netcdf.putAtt(ncid,id_int_time,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_int_time,'VAR_VALIS_MAX','15.0f');
netcdf.putAtt(ncid,id_int_time,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_int_time,'long_name','chirp_integration_time');
netcdf.putAtt(ncid,id_int_time,'short_name','int_time');


%%%%%%%% time dependend variables
id_time = netcdf.defVar(ncid,'DATETIME','nc_uint',did_time);
netcdf.putAtt(ncid,id_time,'VAR_NAME','DATETIME');
netcdf.putAtt(ncid,id_time,'VAR_DESCRIPTION',['Mean time of the measurement '...
                                              'sequence used for each measurement; '...
                                              'defined relative to reference '...
                                              'datetime of Jan. 1 2000 at '...
                                              '0:00:00 UT which is equal to 0.00']);
netcdf.putAtt(ncid,id_time,'VAR_NOTES','');
netcdf.putAtt(ncid,id_time,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_time,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_time,'VAR_DATA_TYPE','DOUBLE');
netcdf.putAtt(ncid,id_time,'VAR_UNITS','MJD2K');
netcdf.putAtt(ncid,id_time,'VAR_SI_CONVERSION','0.0;86400.0;s');
netcdf.putAtt(ncid,id_time,'VAR_VALIS_MIN',nanmin(data.time(:)));
netcdf.putAtt(ncid,id_time,'VAR_VALIS_MAX',nanmax(data.time(:)));
netcdf.putAtt(ncid,id_time,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_time,'long_name','time');
netcdf.putAtt(ncid,id_time,'short_name','time');
netcdf.putAtt(ncid,id_time,'units','seconds since 2001.01.01. 00:00:00.');
if isfield(data, 'totsampchangelabel' )
    netcdf.putAtt(ncid,id_time, 'quality_flag', 'Dublicate time stamps found in lv0-file, the first occurrence of the dublicate time is removed')
end

% id_sampleTms = netcdf.defVar(ncid,'sample_tms','nc_int',did_time);
% netcdf.putAtt(ncid,id_sampleTms,'GEOMS_name','DATETIME.milliseconds');
% netcdf.putAtt(ncid,id_sampleTms,'long_name','Milliseconds of sample');
% netcdf.putAtt(ncid,id_sampleTms,'units','mu s');
% netcdf.putAtt(ncid,id_sampleTms,'comment','To get the correct time the variable sample_tms must be added: time = time + sample_tms.');

id_RR = netcdf.defVar(ncid,'RAIN.RATE.SURFACE','nc_float',did_time);
netcdf.putAtt(ncid,id_RR,'VAR_NAME','RAIN.RATE.SURFACE');
netcdf.putAtt(ncid,id_RR,'VAR_DESCRIPTION',['Rain rate measured by the ' ...
                                            'meteo-station attached to ' ...
                                            'the radar']);
netcdf.putAtt(ncid,id_RR,'VAR_NOTES','');
netcdf.putAtt(ncid,id_RR,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_RR,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_RR,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_RR,'VAR_UNITS','mm h-1');
netcdf.putAtt(ncid,id_RR,'VAR_SI_CONVERSION','0.0;2.777778E-7;m s-1');
netcdf.putAtt(ncid,id_RR,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_RR,'VAR_VALIS_MAX','100.0f');
netcdf.putAtt(ncid,id_RR,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_RR,'long_name','rainfall_rate');
netcdf.putAtt(ncid,id_RR,'short_name','rr');
netcdf.putAtt(ncid,id_RR,'source','Vaisala weather station WXT520 or WXT530');

id_rh = netcdf.defVar(ncid,'HUMIDITY.RELATIVE.SURFACE','nc_float',did_time);
netcdf.putAtt(ncid,id_rh,'VAR_NAME','HUMIDITY.RELATIVE.SURFACE');
netcdf.putAtt(ncid,id_rh,'VAR_DESCRIPTION',['Relative humidity measured by '...
                                            'the meteo-station attached to '...
                                            'the radar']);
netcdf.putAtt(ncid,id_rh,'VAR_NOTES','');
netcdf.putAtt(ncid,id_rh,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_rh,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_rh,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_rh,'VAR_UNITS','%');
netcdf.putAtt(ncid,id_rh,'VAR_SI_CONVERSION','0.0;0.01;1');
netcdf.putAtt(ncid,id_rh,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_rh,'VAR_VALIS_MAX','100.0f');
netcdf.putAtt(ncid,id_rh,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_rh,'long_name','relative_humidity');
netcdf.putAtt(ncid,id_rh,'short_name','rh');
netcdf.putAtt(ncid,id_rh,'source','Vaisala weather station WXT520 or WXT530');

id_T_env = netcdf.defVar(ncid,'SURFACE.TEMPERATURE','nc_float',did_time);
netcdf.putAtt(ncid,id_T_env,'VAR_NAME','SURFACE.TEMPERATURE');
netcdf.putAtt(ncid,id_T_env,'VAR_DESCRIPTION',['Temperature of the environment '...
                                               'measured by the meteo-station '...
                                               'attached to the radar']);
netcdf.putAtt(ncid,id_T_env,'VAR_NOTES','');
netcdf.putAtt(ncid,id_T_env,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_T_env,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_T_env,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_T_env,'VAR_UNITS','K');
netcdf.putAtt(ncid,id_T_env,'VAR_SI_CONVERSION','0.0;1.0;K');
netcdf.putAtt(ncid,id_T_env,'VAR_VALIS_MIN','220.0f');
netcdf.putAtt(ncid,id_T_env,'VAR_VALIS_MAX','330.0f');
netcdf.putAtt(ncid,id_T_env,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_T_env,'long_name','surface_temperature');
netcdf.putAtt(ncid,id_T_env,'short_name','ta');
netcdf.putAtt(ncid,id_T_env,'source','Vaisala weather station WXT520 or WXT530');
netcdf.putAtt(ncid,id_T_env,'comment',['Air temperature is the bulk temperature '...
                                       'of the air, not the surface (skin) temperature.']);

id_pres = netcdf.defVar(ncid,'SURFACE.PRESSURE','nc_float',did_time);
netcdf.putAtt(ncid,id_pres,'VAR_NAME','SURFACE.PRESSURE');
netcdf.putAtt(ncid,id_pres,'VAR_DESCRIPTION',['Surface air pressure '...
                                               'measured by the meteo-station '...
                                               'attached to the radar']);
netcdf.putAtt(ncid,id_pres,'VAR_NOTES','');
netcdf.putAtt(ncid,id_pres,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_pres,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_pres,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_pres,'VAR_UNITS','hPa');
netcdf.putAtt(ncid,id_pres,'VAR_SI_CONVERSION','0.0;1.0E2;kg m-1 s-2');
netcdf.putAtt(ncid,id_pres,'VAR_VALIS_MIN','940.0f');
netcdf.putAtt(ncid,id_pres,'VAR_VALIS_MAX','1050.0f');
netcdf.putAtt(ncid,id_pres,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_pres,'long_name','surface_air_pressure');
netcdf.putAtt(ncid,id_pres,'short_name','pa');
netcdf.putAtt(ncid,id_pres,'source','Vaisala weather station WXT520 or WXT530');
netcdf.putAtt(ncid,id_pres,'comment',['The surface called "surface" means the '...
                                      'lower boundary of the atmosphere. Air '...
                                      'pressure is the force per unit area '...
                                      'which would be exerted when the moving '...
                                      'gas molecules of which the air is '...
                                      'composed strike a theoretical surface '...
                                      'of any orientation.']);

id_ff = netcdf.defVar(ncid,'WIND.SPEED.SURFACE','nc_float',did_time);
netcdf.putAtt(ncid,id_ff,'VAR_NAME','WIND.SPEED.SURFACE');
netcdf.putAtt(ncid,id_ff,'VAR_DESCRIPTION',['Wind speed measured by the '...
                                            'meteo-station attached to the '...
                                            'radar']);
netcdf.putAtt(ncid,id_ff,'VAR_NOTES','');
netcdf.putAtt(ncid,id_ff,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_ff,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_ff,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_ff,'VAR_UNITS','m s-1');
netcdf.putAtt(ncid,id_ff,'VAR_SI_CONVERSION','0.0;1.0;m s-1');
netcdf.putAtt(ncid,id_ff,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_ff,'VAR_VALIS_MAX','100.0f');
netcdf.putAtt(ncid,id_ff,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_ff,'long_name','wind_speed');
netcdf.putAtt(ncid,id_ff,'short_name','wspeed');
netcdf.putAtt(ncid,id_ff,'source','Vaisala weather station WXT520 or WXT530');
netcdf.putAtt(ncid,id_ff,'comment',['Speed is the magnitude of velocity. Wind '...
                                    'is defined as a two-dimensional (horizontal) '...
                                    'air velocity vector, with no vertical '...
                                    'component. (Vertical motion in the atmosphere '...
                                    'has the standard name upward_air_velocity.) '...
                                    'The wind speed is the magnitude of the wind velocity.']);

id_fff = netcdf.defVar(ncid,'WIND.DIRECTION.SURFACE','nc_float',did_time);
netcdf.putAtt(ncid,id_fff,'VAR_NAME','WIND.DIRECTION.SURFACE');
netcdf.putAtt(ncid,id_fff,'VAR_DESCRIPTION',['Wind direction measured by the '...
                                            'meteo-station attached to the '...
                                            'radar']);
netcdf.putAtt(ncid,id_fff,'VAR_NOTES','');
netcdf.putAtt(ncid,id_fff,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_fff,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_fff,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_fff,'VAR_UNITS','deg');
netcdf.putAtt(ncid,id_fff,'VAR_SI_CONVERSION','0.0;1.74533E-2;rad');
netcdf.putAtt(ncid,id_fff,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_fff,'VAR_VALIS_MAX','360.0f');
netcdf.putAtt(ncid,id_fff,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_fff,'long_name','wind_from_direction');
netcdf.putAtt(ncid,id_fff,'short_name','wdir');
netcdf.putAtt(ncid,id_fff,'source','Vaisala weather station WXT520 or WXT530');
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
              
id_RR_source = netcdf.defVar(ncid,'WIND.DIRECTION_SOURCE.SURFACE','nc_float',did_scalar);
netcdf.putAtt(ncid,id_RR_source,'VAR_NAME','WIND.DIRECTION.SURFACE_SOURCE');
netcdf.putAtt(ncid,id_RR_source,'VAR_DESCRIPTION',['Source wind direction measurement'...
                                                   '(Vaisala weather station WXT520 or WXT530)']);
netcdf.putAtt(ncid,id_RR_source,'VAR_NOTES','');
netcdf.putAtt(ncid,id_RR_source,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_RR_source,'VAR_DEPEND','CONSTANT');
netcdf.putAtt(ncid,id_RR_source,'VAR_DATA_TYPE','STRING');
netcdf.putAtt(ncid,id_RR_source,'VAR_UNITS','');
netcdf.putAtt(ncid,id_RR_source,'VAR_SI_CONVERSION','');
netcdf.putAtt(ncid,id_RR_source,'VAR_VALIS_MIN','');
netcdf.putAtt(ncid,id_RR_source,'VAR_VALIS_MAX','');
netcdf.putAtt(ncid,id_RR_source,'VAR_FILL_VALUE','');

id_rh_source = netcdf.defVar(ncid,'HUMIDITY.RELATIVE.SURFACE_SOURCE','nc_float',did_scalar);
netcdf.putAtt(ncid,id_rh_source,'VAR_NAME','HUMIDITY.RELATIVE_SOURCE');
netcdf.putAtt(ncid,id_rh_source,'VAR_DESCRIPTION',['Source relative humidity measurement'...
                                                   '(Vaisala weather station WXT520 or WXT530)']);
netcdf.putAtt(ncid,id_rh_source,'VAR_NOTES','');
netcdf.putAtt(ncid,id_rh_source,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_rh_source,'VAR_DEPEND','CONSTANT');
netcdf.putAtt(ncid,id_rh_source,'VAR_DATA_TYPE','STRING');
netcdf.putAtt(ncid,id_rh_source,'VAR_UNITS','');
netcdf.putAtt(ncid,id_rh_source,'VAR_SI_CONVERSION','');
netcdf.putAtt(ncid,id_rh_source,'VAR_VALIS_MIN','');
netcdf.putAtt(ncid,id_rh_source,'VAR_VALIS_MAX','');
netcdf.putAtt(ncid,id_rh_source,'VAR_FILL_VALUE','');

id_T_env_source = netcdf.defVar(ncid,'SURFACE.TEMPERATURE_SOURCE','nc_float',did_scalar);
netcdf.putAtt(ncid,id_T_env_source,'VAR_NAME','SURFACE.TEMPERATURE_SOURCE');
netcdf.putAtt(ncid,id_T_env_source,'VAR_DESCRIPTION',['Source surface temperature measurement'...
                                                      '(Vaisala weather station WXT520 or WXT530)']);
netcdf.putAtt(ncid,id_T_env_source,'VAR_NOTES','');
netcdf.putAtt(ncid,id_T_env_source,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_T_env_source,'VAR_DEPEND','CONSTANT');
netcdf.putAtt(ncid,id_T_env_source,'VAR_DATA_TYPE','STRING');
netcdf.putAtt(ncid,id_T_env_source,'VAR_UNITS','');
netcdf.putAtt(ncid,id_T_env_source,'VAR_SI_CONVERSION','');
netcdf.putAtt(ncid,id_T_env_source,'VAR_VALIS_MIN','');
netcdf.putAtt(ncid,id_T_env_source,'VAR_VALIS_MAX','');
netcdf.putAtt(ncid,id_T_env_source,'VAR_FILL_VALUE','');

id_pres_source = netcdf.defVar(ncid,'SURFACE.PRESSURE_SOURCE','nc_float',did_scalar);
netcdf.putAtt(ncid,id_pres_source,'VAR_NAME','SURFACE.PRESSURE_SOURCE');
netcdf.putAtt(ncid,id_pres_source,'VAR_DESCRIPTION',['Source surface air pressure measurement'...
                                                     '(Vaisala weather station WXT520 or WXT530)']);
netcdf.putAtt(ncid,id_pres_source,'VAR_NOTES','');
netcdf.putAtt(ncid,id_pres_source,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_pres_source,'VAR_DEPEND','CONSTANT');
netcdf.putAtt(ncid,id_pres_source,'VAR_DATA_TYPE','STRING');
netcdf.putAtt(ncid,id_pres_source,'VAR_UNITS','');
netcdf.putAtt(ncid,id_pres_source,'VAR_SI_CONVERSION','');
netcdf.putAtt(ncid,id_pres_source,'VAR_VALIS_MIN','');
netcdf.putAtt(ncid,id_pres_source,'VAR_VALIS_MAX','');
netcdf.putAtt(ncid,id_pres_source,'VAR_FILL_VALUE','');

id_ff_source = netcdf.defVar(ncid,'WIND.SPEED.SURFACE_SOURCE','nc_float',did_scalar);
netcdf.putAtt(ncid,id_ff_source,'VAR_NAME','WIND.SPEED.SURFACE_SOURCE');
netcdf.putAtt(ncid,id_ff_source,'VAR_DESCRIPTION',['Source wind speed measurement '...
                                                     '(Vaisala weather station WXT520 or WXT530)']);
netcdf.putAtt(ncid,id_ff_source,'VAR_NOTES','');
netcdf.putAtt(ncid,id_ff_source,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_ff_source,'VAR_DEPEND','CONSTANT');
netcdf.putAtt(ncid,id_ff_source,'VAR_DATA_TYPE','STRING');
netcdf.putAtt(ncid,id_ff_source,'VAR_UNITS','');
netcdf.putAtt(ncid,id_ff_source,'VAR_SI_CONVERSION','');
netcdf.putAtt(ncid,id_ff_source,'VAR_VALIS_MIN','');
netcdf.putAtt(ncid,id_ff_source,'VAR_VALIS_MAX','');
netcdf.putAtt(ncid,id_ff_source,'VAR_FILL_VALUE','');


id_fff_source = netcdf.defVar(ncid,'WIND.DIRECTION.SURFACE_SORCE','nc_float',did_scalar);
netcdf.putAtt(ncid,id_fff_source,'VAR_NAME','WIND.DIRECTION.SURFACE_SORCE');
netcdf.putAtt(ncid,id_fff_source,'VAR_DESCRIPTION',['Source wind direction measurement '...
                                                    '(Vaisala weather station WXT520 or WXT530)']);
netcdf.putAtt(ncid,id_fff_source,'VAR_NOTES','');
netcdf.putAtt(ncid,id_fff_source,'VAR_SIZE','1');
netcdf.putAtt(ncid,id_fff_source,'VAR_DEPEND','CONSTANT');
netcdf.putAtt(ncid,id_fff_source,'VAR_DATA_TYPE','STRING');
netcdf.putAtt(ncid,id_fff_source,'VAR_UNITS','');
netcdf.putAtt(ncid,id_fff_source,'VAR_SI_CONVERSION','');
netcdf.putAtt(ncid,id_fff_source,'VAR_VALIS_MIN','');
netcdf.putAtt(ncid,id_fff_source,'VAR_VALIS_MAX','');
netcdf.putAtt(ncid,id_fff_source,'VAR_FILL_VALUE','');
     

id_Tb = netcdf.defVar(ncid,'TEMPERATURE.BRIGHTNESS','nc_float',did_time);
netcdf.putAtt(ncid,id_Tb,'VAR_NAME','TEMPERATURE.BRIGHTNESS');
netcdf.putAtt(ncid,id_Tb,'VAR_DESCRIPTION',['Brightness temperature '...
                                            'at 89-GHz detected by '...
                                            'the passive channel of '...
                                            'the receiver antenna '...
                                            'of the radar']);                                       
netcdf.putAtt(ncid,id_Tb,'VAR_NOTES',['The brightness temperature of a body '...
                                      'is the temperature of a black body '...
                                      'which radiates the same power per '...
                                      'unit solid angle per unit area. Brightness '...
                                      'Temperature measurements at 89-GHz '...
                                      'by the radar']);
netcdf.putAtt(ncid,id_Tb,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_Tb,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_Tb,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_Tb,'VAR_UNITS','K');
netcdf.putAtt(ncid,id_Tb,'VAR_SI_CONVERSION','0.0;1.0;K');
netcdf.putAtt(ncid,id_Tb,'VAR_VALIS_MIN','60.0');
netcdf.putAtt(ncid,id_Tb,'VAR_VALIS_MAX','120.0');
netcdf.putAtt(ncid,id_Tb,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_Tb,'long_name','brightness_temperature');
netcdf.putAtt(ncid,id_Tb,'short_name','tb');


id_lwp = netcdf.defVar(ncid,'LIQUID.WATER.PATH','nc_float',did_time);
netcdf.putAtt(ncid,id_lwp,'VAR_NAME','LIQUID.WATER.PATH');
netcdf.putAtt(ncid,id_lwp,'VAR_DESCRIPTION',['Retrieval based on Nural '...
                                             'NetWork developed by RPG '...
                                             'and trained on Radiosondes '...
                                             'or NWP reanalysis data']);
netcdf.putAtt(ncid,id_lwp,'VAR_NOTES',['A measure of the weight of the liquid '...
                                     'water droplets in the atmosphere above '...
                                     'a unit surface area on the earth. The '...
                                     'liquid  water path is calculated from '...
                                     'the tb measurement of the 89-GHz chanal. '...
                                     'The retrieval was developed by RPG and '...
                                     'is based on a neural network approach']);
netcdf.putAtt(ncid,id_lwp,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_lwp,'VAR_DEPEND','REAL');
netcdf.putAtt(ncid,id_lwp,'VAR_DATA_TYPE','DATETIME');
netcdf.putAtt(ncid,id_lwp,'VAR_UNITS','g m-2');
netcdf.putAtt(ncid,id_lwp,'VAR_SI_CONVERSION','0.0;1.0E-3;kg m-2');
netcdf.putAtt(ncid,id_lwp,'VAR_VALIS_MIN','0.0');
netcdf.putAtt(ncid,id_lwp,'VAR_VALIS_MAX','1000.0');
netcdf.putAtt(ncid,id_lwp,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_lwp,'long_name','liquid_water_path');
netcdf.putAtt(ncid,id_lwp,'short_name','lwp');

id_TransPow = netcdf.defVar(ncid,'TRANSMITTER.POWER','nc_float',did_time);
netcdf.putAtt(ncid,id_TransPow,'VAR_NAME','TRANSMITTER.POWER');
netcdf.putAtt(ncid,id_TransPow,'VAR_DESCRIPTION','Radar Transmitter Power');
netcdf.putAtt(ncid,id_TransPow,'VAR_NOTES',['The transmitted power of the '...
                                            'radar instrument is part of '...
                                            'the house keeping data and used '...
                                            'to monitor its the performance']);
netcdf.putAtt(ncid,id_TransPow,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_TransPow,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_TransPow,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_TransPow,'VAR_UNITS','W');
netcdf.putAtt(ncid,id_TransPow,'VAR_SI_CONVERSION','0.0;1.0;m2 kg s-3');
netcdf.putAtt(ncid,id_TransPow,'VAR_VALIS_MIN','0.0');
netcdf.putAtt(ncid,id_TransPow,'VAR_VALIS_MAX','10.0');
netcdf.putAtt(ncid,id_TransPow,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_TransPow,'long_name','platform_transmitter_power');
netcdf.putAtt(ncid,id_TransPow,'short_name','p_trans');

                                 
id_T_trans = netcdf.defVar(ncid,'TRANSMITTER.TEMPERATURE','nc_float',did_time);
netcdf.putAtt(ncid,id_T_trans,'VAR_NAME','TRANSMITTER.TEMPERATURE');
netcdf.putAtt(ncid,id_T_trans,'VAR_DESCRIPTION','Radar Transmitter Temperature');
netcdf.putAtt(ncid,id_T_trans,'VAR_NOTES',['The transmitter temperature of '...
                                           'the radar instrument is part of '...
                                           'the house keeping data and used '...
                                           'to monitor its the performance. '...
                                           'For a stable performance 309K < '...
                                           't_trans < 313K.']);
netcdf.putAtt(ncid,id_T_trans,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_T_trans,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_T_trans,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_T_trans,'VAR_UNITS','K');
netcdf.putAtt(ncid,id_T_trans,'VAR_SI_CONVERSION','0.0;1.0;K');
netcdf.putAtt(ncid,id_T_trans,'VAR_VALIS_MIN','280.0');
netcdf.putAtt(ncid,id_T_trans,'VAR_VALIS_MAX','330.0');
netcdf.putAtt(ncid,id_T_trans,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_T_trans,'long_name','platform_transmitter_temperature');
netcdf.putAtt(ncid,id_T_trans,'short_name','t_trans');


id_T_rec = netcdf.defVar(ncid,'RECEIVER.TEMPERATURE','nc_float',did_time);
netcdf.putAtt(ncid,id_T_rec,'VAR_NAME','RECEIVER.TEMPERATURE');
netcdf.putAtt(ncid,id_T_rec,'VAR_DESCRIPTION','Radar Receiver Temperature');
netcdf.putAtt(ncid,id_T_rec,'VAR_NOTES',['The receiver temperature of the radar '...
                                         'instrument is part of the house '...
                                         'keeping data and used to monitor '...
                                         'its the performance. For a stable '...
                                         'performance 306K < t_rec < 312K.']);
netcdf.putAtt(ncid,id_T_rec,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_T_rec,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_T_rec,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_T_rec,'VAR_UNITS','K');
netcdf.putAtt(ncid,id_T_rec,'VAR_SI_CONVERSION','0.0;1.0;K');
netcdf.putAtt(ncid,id_T_rec,'VAR_VALIS_MIN','280.0');
netcdf.putAtt(ncid,id_T_rec,'VAR_VALIS_MAX','330.0');
netcdf.putAtt(ncid,id_T_rec,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_T_rec,'long_name','platform_receiver_temperature');
netcdf.putAtt(ncid,id_T_rec,'short_name','t_rec');
                                 

id_T_pc = netcdf.defVar(ncid,'COMPUTER.TEMPERATURE','nc_float',did_time);
netcdf.putAtt(ncid,id_T_rec,'VAR_NAME','COMPUTER.TEMPERATURE');
netcdf.putAtt(ncid,id_T_rec,'VAR_DESCRIPTION','Temperature of the radar controlling PC');
netcdf.putAtt(ncid,id_T_rec,'VAR_NOTES',['The radar controlling PC temperature '...
                                         'is part of the house keeping data '...
                                         'and used to monitor its the '...
                                         'performance. For a stable performance '...
                                         't_pc < 323K.']);
netcdf.putAtt(ncid,id_T_rec,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_T_rec,'VAR_DEPEND','DATETIME');
netcdf.putAtt(ncid,id_T_rec,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_T_rec,'VAR_UNITS','K');
netcdf.putAtt(ncid,id_T_rec,'VAR_SI_CONVERSION','0.0;1.0;K');
netcdf.putAtt(ncid,id_T_rec,'VAR_VALIS_MIN','270.0');
netcdf.putAtt(ncid,id_T_rec,'VAR_VALIS_MAX','350.0');
netcdf.putAtt(ncid,id_T_rec,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_T_pc,'long_name','platform_computer_temperature');
netcdf.putAtt(ncid,id_T_pc,'short_name','t_pc');

                                 
%%%%%%%% multi-D variables

id_Ze = netcdf.defVar(ncid,'RADAR.REFLECTIVITY.FACTOR','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_Ze,'VAR_NAME','RADAR.REFLECTIVITY.FACTOR');
netcdf.putAtt(ncid,id_Ze,'VAR_DESCRIPTION',['Equivalent radar reflectivity ' ...
                                            'factor at vertical polarisation']);
netcdf.putAtt(ncid,id_Ze,'VAR_NOTES',['Equivalent reflectivity factor is the '...
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
netcdf.putAtt(ncid,id_Ze,'VAR_SIZE','2;RANGE');
netcdf.putAtt(ncid,id_Ze,'VAR_DEPEND','DATETIME;RANGE');
netcdf.putAtt(ncid,id_Ze,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_Ze,'VAR_UNITS','mm6 m-3');
netcdf.putAtt(ncid,id_Ze,'VAR_SI_CONVERSION','0.0;1.0E-18;m3');
netcdf.putAtt(ncid,id_Ze,'VAR_VALIS_MIN','1.e-10');
netcdf.putAtt(ncid,id_Ze,'VAR_VALIS_MAX','10.f');
netcdf.putAtt(ncid,id_Ze,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_Ze,'long_name','equivalent_reflectivity_factor');
netcdf.putAtt(ncid,id_Ze,'short_name','ze');
if isfield(data, 'Ze_label') % Ze corrected, adding note
    netcdf.putAtt(ncid,id_Ze,'comment2',data.Ze_label);
    netcdf.putAtt(ncid,id_Ze,'corretion_dB',data.Ze_corr);
end

id_vm = netcdf.defVar(ncid,'DOPPLER.VELOCITY_MEAN','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_vm,'VAR_NAME','DOPPLER.VELOCITY_MEAN');
netcdf.putAtt(ncid,id_vm,'VAR_DESCRIPTION',['Mean Doppler velocity at ' ...
                                            'vertical polarisation']);
netcdf.putAtt(ncid,id_vm,'VAR_NOTES',['The Doppler velocity is the radial component of '...
                                    'the velocity vector of a scattering object as '...
                                    'observed by a remote sensor, such as a Doppler radar. '...
                                    'NOTE: Negative velocities indicate particles '...
                                    'motion towards the radar']);
netcdf.putAtt(ncid,id_vm,'VAR_SIZE','2;RANGE');
netcdf.putAtt(ncid,id_vm,'VAR_DEPEND','DATETIME;RANGE');
netcdf.putAtt(ncid,id_vm,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_vm,'VAR_UNITS','m s-1');
netcdf.putAtt(ncid,id_vm,'VAR_SI_CONVERSION','0.0;1.0;m s-1');
netcdf.putAtt(ncid,id_vm,'VAR_VALIS_MIN','0.0');
netcdf.putAtt(ncid,id_vm,'VAR_VALIS_MAX','15.0');
netcdf.putAtt(ncid,id_vm,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_vm,'long_name','mean_doppler_velocity');
netcdf.putAtt(ncid,id_vm,'short_name','vm');

id_sigma = netcdf.defVar(ncid,'DOPPLER.SPECTRUM_WIDTH','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_sigma,'VAR_NAME','DOPPLER.SPECTRUM_WIDTH');
netcdf.putAtt(ncid,id_sigma,'VAR_DESCRIPTION',['Doppler spectrum width at' ... 
                                               'vertical polarisation']);
netcdf.putAtt(ncid,id_sigma,'VAR_NOTES','');
netcdf.putAtt(ncid,id_sigma,'VAR_SIZE','2;RANGE');
netcdf.putAtt(ncid,id_sigma,'VAR_DEPEND','DATETIME;RANGE');
netcdf.putAtt(ncid,id_sigma,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_sigma,'VAR_UNITS','m s-1');
netcdf.putAtt(ncid,id_sigma,'VAR_SI_CONVERSION','0.0;1.0;m s-1');
netcdf.putAtt(ncid,id_sigma,'VAR_VALIS_MIN','0.0');
netcdf.putAtt(ncid,id_sigma,'VAR_VALIS_MAX','5.0');
netcdf.putAtt(ncid,id_sigma,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_sigma,'long_name','doppler_spectrum_width');
netcdf.putAtt(ncid,id_sigma,'short_name','sw');

id_skew = netcdf.defVar(ncid,'DOPPLER.SPECTRUM_SKEWNESS','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_skew,'VAR_NAME','DOPPLER.SPECTRUM_SKEWNESS');
netcdf.putAtt(ncid,id_skew,'VAR_DESCRIPTION',['Doppler spectrum skewness at' ...
                                              'vertical polarisation']);
netcdf.putAtt(ncid,id_skew,'VAR_NOTES','');
netcdf.putAtt(ncid,id_skew,'VAR_SIZE','2;RANGE');
netcdf.putAtt(ncid,id_skew,'VAR_DEPEND','DATETIME;RANGE');
netcdf.putAtt(ncid,id_skew,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_skew,'VAR_UNITS','1');
netcdf.putAtt(ncid,id_skew,'VAR_SI_CONVERSION','0.0;1.0;1');
netcdf.putAtt(ncid,id_skew,'VAR_VALIS_MIN','-100.0');
netcdf.putAtt(ncid,id_skew,'VAR_VALIS_MAX','100.0');
netcdf.putAtt(ncid,id_skew,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_skew,'long_name','doppler_spectrum_skewnes');
netcdf.putAtt(ncid,id_skew,'short_name','skew');

%Included by Bravo-Aranda, J.A. JABA
if data.DualPol > 0
    id_Ze_hv = netcdf.defVar(ncid,'ze_hv','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_Ze_hv,'GEOMS_name','RADAR.REFLECTIVITY.FACTOR_HV');
    netcdf.putAtt(ncid,id_Ze_hv,'long_name','equivalent_reflectivity_factor');
    netcdf.putAtt(ncid,id_Ze_hv,'rpg_name','Ze_hv');
    netcdf.putAtt(ncid,id_Ze_hv,'standard_name','Equivalent radar reflectivity factor of the cross polarized channel');
    netcdf.putAtt(ncid,id_Ze_hv,'units','mm6 m-3');
    netcdf.putAtt(ncid,id_Ze_hv,'valid_range',[nanmin(data.Ze_hv(:)), nanmax(data.Ze_hv(:))]);
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
    
    id_ldr = netcdf.defVar(ncid,'LINEAR.DEPOLARIZATION.RATIO','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_ldr,'VAR_NAME','LINEAR.DEPOLARIZATION.RATIO');
    netcdf.putAtt(ncid,id_ldr,'VAR_DESCRIPTION',['Doppler spectrum skewness ' ...
                                                 'at vertical polarisation']);
    netcdf.putAtt(ncid,id_ldr,'VAR_NOTES',['The ratio of the power received in '...
                                         'the orthogonal, or cross-polarized, '...
                                         'channel to that received in the '...
                                         'transmission, or copolarized, channel'...
                                         'of a dual- channel radar, when a '...
                                         'linearly polarized signal is transmitted. '...
                                         'L_dr = Ze_hv/Ze_vv']);
    netcdf.putAtt(ncid,id_ldr,'VAR_SIZE','2;RANGE');
    netcdf.putAtt(ncid,id_ldr,'VAR_DEPEND','DATETIME;RANGE');
    netcdf.putAtt(ncid,id_ldr,'VAR_DATA_TYPE','REAL');
    netcdf.putAtt(ncid,id_ldr,'VAR_UNITS','mm6 m-3');
    netcdf.putAtt(ncid,id_ldr,'VAR_SI_CONVERSION','0.0;1.0E-18;m3');
    netcdf.putAtt(ncid,id_ldr,'VAR_VALIS_MIN','1.e-09f');
    netcdf.putAtt(ncid,id_ldr,'VAR_VALIS_MAX','100000.f');
    netcdf.putAtt(ncid,id_ldr,'VAR_FILL_VALUE',fill_value);    
    netcdf.putAtt(ncid,id_ldr,'long_name','linear_depolarisation_ratio');
    netcdf.putAtt(ncid,id_ldr,'short_name','ldr');
end

if data.DualPol == 2
    id_xcorr = netcdf.defVar(ncid,'rho_hv','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_xcorr,'GEOMS_name','CORRELATION.COEFFICIENT.');
    netcdf.putAtt(ncid,id_xcorr,'long_name','correlation_coefficient');
    netcdf.putAtt(ncid,id_xcorr,'rpg_name','CorrCoeff');    
    netcdf.putAtt(ncid,id_xcorr,'standard_name','co-cross-channel correlation coefficient');
    netcdf.putAtt(ncid,id_xcorr,'short_name','rho_hv');        
%     netcdf.putAtt(ncid,id_xcorr,'valid_range',[nanmin(data.xcorr(:)), nanmax(data.xcorr(:))]);
    netcdf.putAtt(ncid,id_xcorr,'fill_value',str_fill_value);
    
    id_difphase = netcdf.defVar(ncid,'phi_dp','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_difphase,'GEOMS_name','DIFFERENTIAL.PHASE');
    netcdf.putAtt(ncid,id_difphase,'long_name','differential_phase');
    netcdf.putAtt(ncid,id_difphase,'rpg_name','difphase');    
    netcdf.putAtt(ncid,id_difphase,'standard_name','co-cross-channel differential phase');   
    netcdf.putAtt(ncid,id_difphase,'short_name','phi_dp');        
    netcdf.putAtt(ncid,id_difphase,'unite','degree');   
%     netcdf.putAtt(ncid,id_difphase,'valid_range',[nanmin(data.difphase(:)), nanmax(data.difphase(:))]); 
    netcdf.putAtt(ncid,id_difphase,'fill_value',str_fill_value); 
end

id_QualFlag = netcdf.defVar(ncid,'FLAG.PROCESSING.QUALITY','nc_byte',[did_range,did_time]);
netcdf.putAtt(ncid,id_QualFlag,'VAR_NAME','FLAG.PROCESSING.QUALITY');
netcdf.putAtt(ncid,id_QualFlag,'VAR_DESCRIPTION',['If 2^0 bit is 1: this range gate ' ...
                                            'is known to have aritifical spikes ', ...
                                            'occurring. If 2^1 bit is 1: ' ...
                                            'aircraft or other known flying ' ...
                                            'non-meteorological object ', ...
                                            'If 2^2 bit is 1: wet-radome']);
netcdf.putAtt(ncid,id_QualFlag,'VAR_NOTES',['This variable contains information '...
                                            'on anything that might impact ' ...
                                            'the quality of the data at each ' ...
                                            'pixel. Must be converted into ' ...
                                            'three bit binary string. ', ...
                                            'If 0, i.e. dec2bin(QualityFlag,3) ' ...
                                            '= 000, none of the included issues ', ...
                                            'were found. The definitions of ' ...
                                            'each bit are given in the definition ' ... 
                                            'attribute.']);
netcdf.putAtt(ncid,id_QualFlag,'VAR_SIZE','2;RANGE');
netcdf.putAtt(ncid,id_QualFlag,'VAR_DEPEND',['DATETIME;',len_range]);
netcdf.putAtt(ncid,id_QualFlag,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_QualFlag,'VAR_UNITS','1');
netcdf.putAtt(ncid,id_QualFlag,'VAR_SI_CONVERSION','0.0;1.0;1');
netcdf.putAtt(ncid,id_QualFlag,'VAR_VALIS_MIN','0.0');
netcdf.putAtt(ncid,id_QualFlag,'VAR_VALIS_MAX','10.0');
netcdf.putAtt(ncid,id_QualFlag,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_QualFlag,'long_name','quality_flag_data_processing');
netcdf.putAtt(ncid,id_QualFlag,'short_name','pro_qf');

id_Aliasmask = netcdf.defVar(ncid,'ALIAS_MASK','nc_byte',[did_range,did_time]);
netcdf.putAtt(ncid,id_Aliasmask,'VAR_NAME','ALIAS_MASK');
netcdf.putAtt(ncid,id_Aliasmask,'VAR_DESCRIPTION',['Mask indicating in which '...
                                                   'bin dealiaising was applyed. '...
                                                   'This indication is related to '...
                                                   'the processing applyed to ' ...
                                                   'the data.']);
netcdf.putAtt(ncid,id_Aliasmask,'VAR_NOTES',['AliasMask indicates the bins ' ...
                                             'where aliasing is detected. '...
                                             '0 = no aliasing; 1 = aliasing ' ...
                                             'occurs. Only applicable if '...
                                             'variable alias_flag is 2.']);
netcdf.putAtt(ncid,id_Aliasmask,'VAR_SIZE',['2;',len_range]);
netcdf.putAtt(ncid,id_Aliasmask,'VAR_DEPEND','DATETIME;RANGE');
netcdf.putAtt(ncid,id_Aliasmask,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_Aliasmask,'VAR_UNITS','1');
netcdf.putAtt(ncid,id_Aliasmask,'VAR_SI_CONVERSION','0.0;1.0;1');
netcdf.putAtt(ncid,id_Aliasmask,'VAR_VALIS_MIN','0.0');
netcdf.putAtt(ncid,id_Aliasmask,'VAR_VALIS_MAX','10.0');
netcdf.putAtt(ncid,id_Aliasmask,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_Aliasmask,'long_name','alaising_mask');


id_incel = netcdf.defVar(ncid,'ANGLE_PLATFORM.ROLL','nc_float',did_time); 
netcdf.putAtt(ncid,id_incel,'VAR_NAME','ANGLE_PLATFORM.ROLL');
netcdf.putAtt(ncid,id_incel,'VAR_DESCRIPTION',['Roll is a rotation about an ' ...
                                               'axis that is perpendicular ' ...
                                               'to the local vertical axis ' ...
                                               'and is coplanar with the nominal ' ...
                                               'forward motion direction of ' ...
                                               'the platform. Because the radar ' ...
                                               'is fixed to the ground the ' ...
                                               'forward motion axis is not existing. ' ...
                                               'Therefore  the vertical axis ' ...
                                               '(longer dimension of the radar) ' ...
                                               'is set as the moving direction']);
netcdf.putAtt(ncid,id_incel,'VAR_NOTES','');
netcdf.putAtt(ncid,id_incel,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_incel,'VAR_DEPEND','DATETIME;');
netcdf.putAtt(ncid,id_incel,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_incel,'VAR_UNITS','1');
netcdf.putAtt(ncid,id_incel,'VAR_SI_CONVERSION','0.0;1.74533E-2;rad');
netcdf.putAtt(ncid,id_incel,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_incel,'VAR_VALIS_MAX','360.0f');
netcdf.putAtt(ncid,id_incel,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_incel,'long_name','platform_roll_angle');
                            

id_incea = netcdf.defVar(ncid,'ANGLE_PLATFORM.PITCH','nc_float',did_time);
netcdf.putAtt(ncid,id_incea,'VAR_NAME','ANGLE_PLATFORM.PITCH');
netcdf.putAtt(ncid,id_incea,'VAR_DESCRIPTION',['Pitch is a rotation about an ' ...
                                               'axis that is perpendicular ' ...
                                               'to both the local vertical ' ...
                                               'axis and the nominal forward ' ...
                                               'motion direction of the platform. ' ...
                                               'Because the radar is fixed to ' ...
                                               'the ground the forward motion ' ...
                                               'axis is not existing. Therefore ' ...
                                               'the vertical axis (longer dimension ' ...
                                               'of the radar) is set as the ' ...
                                               'moving direction']);
netcdf.putAtt(ncid,id_incea,'VAR_NOTES','');
netcdf.putAtt(ncid,id_incea,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_incea,'VAR_DEPEND','DATETIME;');
netcdf.putAtt(ncid,id_incea,'VAR_DATA_TYPE','REAL');
netcdf.putAtt(ncid,id_incea,'VAR_UNITS','1');
netcdf.putAtt(ncid,id_incea,'VAR_SI_CONVERSION','0.0;1.74533E-2;rad');
netcdf.putAtt(ncid,id_incea,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_incea,'VAR_VALIS_MAX','360.0f');
netcdf.putAtt(ncid,id_incea,'VAR_FILL_VALUE',fill_value);
netcdf.putAtt(ncid,id_incea,'long_name','platform_pitch_angle');


id_QF = netcdf.defVar(ncid,'FLAG.MEASUREMENT.QUALITY','nc_byte',did_time);
netcdf.putAtt(ncid,id_QF,'VAR_NAME','ANGLE_PLATFORM.PITCH');
netcdf.putAtt(ncid,id_QF,'VAR_DESCRIPTION',['bit6 = ADC saturation, ' ...
                                            'bit5 = spectral width too high, ' ...
                                            'bit4 = no transmitter power leveling. ' ...
                                            'bit1 = heater off '...
                                            'bit0 = blower on. '...
                                            'NOTE: convention holds: ' ...
                                            'bit0 = 2^5, bit1 = 2^4,' ...
                                            'bit4 = 2^2, bit5 = 2^1, bit6 = 2^0']);
%                                      'convention holds: bit1 = 2^3, '...
%                                      'bit2 = 2^2, bit3 = 2^1, bit4 = 2^0'
                                 
netcdf.putAtt(ncid,id_QF,'VAR_NOTES',['Quality flag related to the radar '... 
                                      'operation given by radar software ' ...
                                      'from RPG and the hearer and blower status']);
netcdf.putAtt(ncid,id_QF,'VAR_SIZE','2');
netcdf.putAtt(ncid,id_QF,'VAR_DEPEND','DATETIME;');
netcdf.putAtt(ncid,id_QF,'VAR_DATA_TYPE','BINARY');
netcdf.putAtt(ncid,id_QF,'VAR_UNITS','1');
netcdf.putAtt(ncid,id_QF,'VAR_SI_CONVERSION','0.0;1.0;1');
netcdf.putAtt(ncid,id_QF,'VAR_VALIS_MIN','0.0f');
netcdf.putAtt(ncid,id_QF,'VAR_VALIS_MAX','128.0f');
netcdf.putAtt(ncid,id_QF,'VAR_FILL_VALUE',fill_value);

id_status = netcdf.defVar(ncid,'blower_status','nc_byte',did_time);
netcdf.putAtt(ncid,id_status,'long_name','platform_blower_status');
netcdf.putAtt(ncid,id_status,'short_name','blowe_status');
netcdf.putAtt(ncid,id_status,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_status,'comment',['The quality flag shows the following '...
                                        'status: 0/1 = heater on/off; '...
                                        '0/10 = blower on/off. The parameter '...
                                        'is recorded to check and monitor '...
                                        'the radar performance']);           

% id_QF = netcdf.defVar(ncid,'rad_qf','nc_byte',did_time);
% netcdf.putAtt(ncid,id_QF,'GEOMS_name','FLAG.MEASUREMENT.QUALITY');
% netcdf.putAtt(ncid,id_QF,'long_name','quality_falg_platform_measurement');
% netcdf.putAtt(ncid,id_QF,'rpg_name','QF');
% netcdf.putAtt(ncid,id_QF,'standard_name',['Quality flag related to the radar operation '... 
%                                       'given by radar software from RPG']);
% netcdf.putAtt(ncid,id_QF,'short_name','rad_qf');                                  
% netcdf.putAtt(ncid,id_QF,'comment', ['To get the bit entries, one has to'...
%                                      'convert the integer into a 4 bit binary. '...
%                                      'bit4 = ADC saturation, bit3 = spectral '...
%                                      'width too high, bit2 = no transmitter '...
%                                      'power leveling. Note that in the above '...
%                                      'convention holds: bit1 = 2^3, '...
%                                      'bit2 = 2^2, bit3 = 2^1, bit4 = 2^0'])
                                    
                                    
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
netcdf.putVar(ncid, id_lat, 0, data.Lat);
netcdf.putVar(ncid, id_lon, 0, data.Lon);
netcdf.putVar(ncid, id_MSL, 0, data.MSL);
netcdf.putVar(ncid, id_freq, 0, 299792458/(data.freq * 1e9) ); 
netcdf.putVar(ncid, id_wl, 0, data.freq * 1e9);
netcdf.putVar(ncid, id_HPBW, 0, data.HPBW);

temporary_data = 1;    %no need for terminator on longest lines
netcdf.putVar(ncid, id_RR_source,    0, temporary_data);
netcdf.putVar(ncid, id_rh_source,    0, temporary_data);
netcdf.putVar(ncid, id_T_env_source, 0, temporary_data);
netcdf.putVar(ncid, id_pres_source,  0, temporary_data);
netcdf.putVar(ncid, id_ff_source,    0, temporary_data);
netcdf.putVar(ncid, id_fff_source,   0, temporary_data);
clear temporary_data ;

% range dependet
netcdf.putVar(ncid,id_range,0,data.n_levels,data.range);

% chrip seq dependent variables
netcdf.putVar(ncid,id_range_offsets,0,data.no_chirp_seq,data.range_offsets);
vm_res = NaN(size(data.range_offsets)) ;
for ii = 1:length(data.range_offsets) 
    vm_res(ii) = (data.velocity(ii,2)-data.velocity(ii,1)) ;
end 
netcdf.putVar(ncid,id_vm_res,0,data.no_chirp_seq,vm_res);
temporary_data = (1./((4.*data.DoppMax.*94e+10)/(3.0e+09))).*double(data.SeqAvg) ;
netcdf.putVar(ncid,id_int_time,0,data.no_chirp_seq,temporary_data);
clear temporary_data
netcdf.putVar(ncid,id_DoppMax,0,data.no_chirp_seq,data.DoppMax);
netcdf.putVar(ncid,id_AntiAlias,0,data.AntiAlias);

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
clear temporary_data ;
data_qf                      = data.QF ;
data_status                  = data.status ;
data_status(data_status==1)  = 16 ;
data_status(data_status==10) = 32 ;
data_status(data_status==11) = 32+16 ;
temporary_data               = int8(data_status)+data_qf ;
%temporary_data               = dec2bin(temporary_data,6) ;
netcdf.putVar(ncid,id_QF,0,data.totsamp,temporary_data) ;
clear data_qf
clear data_status
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
%     clear temporary_data ;
%     temporary_data = data.Ze_hv ;
%     temporary_data(isnan(temporary_data))=fill_value ;
%     netcdf.putVar(ncid,id_Ze_hv,[0,0],[data.n_levels,data.totsamp],temporary_data'); %LP
end

if data.DualPol == 2
    disp('WARNING! skipping writing id_xcorr and difphase into files, variables not available (yet??)')
    % netcdf.putVar(ncid,id_xcorr,[0,0],[data.n_levels,data.totsamp],data.xcorr'); %JABA
    % netcdf.putVar(ncid,id_difphase,[0,0],[data.n_levels,data.totsamp],data.difphase'); %JABA
end

netcdf.close(ncid);

end