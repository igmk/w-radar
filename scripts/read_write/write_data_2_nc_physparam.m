function write_data_2_nc_physparam(data, outfile, config)
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
did_no_seq = netcdf.defDim(ncid,'chirp',data.no_chirp_seq);
did_scalar = netcdf.defDim(ncid,'scalar',1);

if ne(data.T_altcount,0) && isfield(data,'Tprof')
    did_T_range = netcdf.defDim(ncid,'T_range',data.T_altcount);
end

if ne(data.H_altcount,0) && isfield(data,'Hprof')
    did_H_range = netcdf.defDim(ncid,'H_range',data.H_altcount);
end


%% ######################## add global attributes
glob = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,glob,'FillValue',str_fill_value);
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


% variables for dimensions
id_scal = netcdf.defVar(ncid,'scalar','nc_int',did_scalar);
netcdf.putAtt(ncid,id_scal,'comment','Variable for scalar-dimension');

id_no_seq = netcdf.defVar(ncid,'chirp','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_no_seq,'comment','Variable for chirp-dimension');

id_veldim = netcdf.defVar(ncid,'velocity','nc_int',did_vel);
netcdf.putAtt(ncid,id_veldim,'comment','Variable for velocity-dimension');


%% ################ get variable ids and add attributes

%%%%%%%%%% scalar variables

id_samples = netcdf.defVar(ncid,'samples','nc_int',did_scalar);
netcdf.putAtt(ncid,id_samples,'standard_name','No. samples');

id_levels = netcdf.defVar(ncid,'n_levels','nc_int',did_scalar);
netcdf.putAtt(ncid,id_levels,'standard_name','Number of range levels');

id_no_chirp_seq = netcdf.defVar(ncid,'no_chirp_seq','nc_int',did_scalar);
netcdf.putAtt(ncid,id_no_chirp_seq,'standard_name','Number of chirp sequences');
netcdf.putAtt(ncid,id_no_chirp_seq,'comment',...
    'The radar can be programmed to run different resolution modes for different layers spanning several range gates');



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
                                                                    
if exist('did_T_range', 'var')
    id_T_levels = netcdf.defVar(ncid,'T_levels','nc_float',did_T_range);
    netcdf.putAtt(ncid,id_T_levels,'standard_name','Levels of temperature profile');
end

if exist('did_H_range', 'var')
    id_H_levels = netcdf.defVar(ncid,'H_levels','nc_float',did_H_range);
    netcdf.putAtt(ncid,id_H_levels,'standard_name','Levels of humidity profiles');
end

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


id_dr = netcdf.defVar(ncid,'dz_chirp','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_dr,'GEOMS_name','radar.operation.parameter.dz_chirp');
netcdf.putAtt(ncid,id_dr,'standard_name','Range resolution for chirp sequences');
netcdf.putAtt(ncid,id_dr,'units','m');


id_dv = netcdf.defVar(ncid,'dv_chirp','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_dv,'standard_name','Doppler velocity resolution for chirp sequences');
netcdf.putAtt(ncid,id_dv,'units','m s-1');
netcdf.putAtt(ncid,id_dv,'comment','Calculated from Nyquist velocity and length of Doppler array (DoppRes = 2*nqv/dopp_len).');



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
netcdf.putAtt(ncid,id_rh,'long_name','relative_humidity');
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
netcdf.putAtt(ncid,id_fff,'standard_name','wind direction measured by the meteo-station attached to the radar');
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
                                      'by the passive channel of the receiver '...
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


if isfield(config,'ignoreLWP') % check if field exist
    Flag_noLWP = 1;
    if config.ignoreLWP == 0 % also ignoreLWP = 0 should be taken to indicate to include LWP
        Flag_noLWP = 0;
    end
else
    Flag_noLWP = 0;
end

if ~Flag_noLWP

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
                                         'is based on a neural network approach']);
end

id_sampleTms = netcdf.defVar(ncid,'sample_tms','nc_int',did_time);
netcdf.putAtt(ncid,id_sampleTms,'GEOMS_name','DATETIME.MILLISECONDS');
netcdf.putAtt(ncid,id_sampleTms,'standard_name','Milliseconds of sample');
netcdf.putAtt(ncid,id_sampleTms,'units','mu s');
netcdf.putAtt(ncid,id_sampleTms,'comment','To get the correct time the variable sample_tms must be added: time = time + sample_tms.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% multi-D variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

id_vel = netcdf.defVar(ncid,'Doppler_velocity','nc_float',[did_vel,did_range,did_time]);
netcdf.putAtt(ncid,id_vel,'standard_name','Doppler velocity array for each time-range-bin');
netcdf.putAtt(ncid,id_vel,'units','m/s');
netcdf.putAtt(ncid,id_vel,'comment',...
    ['Doppler velocity array for each bin. Note that if spectra has been dealiased, the ' ...
    'Doppler velocity array is changing from bin to bin']);



%--- radar moments ---

id_Ze = netcdf.defVar(ncid,'ze','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_Ze,'GEOMS_name','RADAR.REFLECTIVITY.FACTOR_VV');
netcdf.putAtt(ncid,id_Ze,'long_name','equivalent_reflectivity_factor');
netcdf.putAtt(ncid,id_Ze,'rpg_name','Ze');
if data.DualPol == 0
    netcdf.putAtt(ncid,id_Ze,'standard_name','Equivalent radar reflectivity factor');
else
    netcdf.putAtt(ncid,id_Ze,'standard_name','Equivalent radar reflectivity factor at vertical polarisation');
end
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
if data.DualPol == 0
    netcdf.putAtt(ncid,id_vm,'standard_name','Mean Doppler velocity');
else
    netcdf.putAtt(ncid,id_vm,'standard_name','Mean Doppler velocity at vertical polarisation');
end
netcdf.putAtt(ncid,id_vm,'short_name','vm');
netcdf.putAtt(ncid,id_vm,'units','m s-1');
netcdf.putAtt(ncid,id_vm,'valid_range',[nanmin(data.vm(:)), nanmax(data.vm(:))]);
netcdf.putAtt(ncid,id_vm,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_vm,'comment',['The Doppler velocity is the radial component of '...
                                    'the velocity vector of a scattering object as '...
                                    'observed by a remote sensor, such as a Doppler radar. '...
                                    'NOTE: Negative velocities indicate particles '...
                                    'motion towards the radar'])

id_sigma = netcdf.defVar(ncid,'sw','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_sigma,'GEOMS_name','DOPPLER.SPECTRUM_WIDTH');
netcdf.putAtt(ncid,id_sigma,'long_name','doppler_spectrum_width');
netcdf.putAtt(ncid,id_sigma,'rpg_name','sigma');
if data.DualPol == 0
    netcdf.putAtt(ncid,id_sigma,'standard_name','Doppler spectrum width');
else
    netcdf.putAtt(ncid,id_sigma,'standard_name','Doppler spectrum width at vertical polarisation');
end
netcdf.putAtt(ncid,id_sigma,'short_name','sw');
netcdf.putAtt(ncid,id_sigma,'units','m s-1');
netcdf.putAtt(ncid,id_sigma,'valid_range',[nanmin(data.sigma(:)), nanmax(data.sigma(:))]);
netcdf.putAtt(ncid,id_sigma,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_sigma,'comment',['Doppler spectrum width of the measured Doppler spectrum'...
                                       'at vertical polarisation.']);

id_skew = netcdf.defVar(ncid,'sk','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_skew,'GEOMS_name','DOPPLER.SPECTRUM_SKEWNESS');
netcdf.putAtt(ncid,id_skew,'long_name','doppler_spectrum_skewnes');
netcdf.putAtt(ncid,id_skew,'rpg_name','skew');
if data.DualPol == 0
    netcdf.putAtt(ncid,id_skew,'standard_name','Doppler spectrum skewness');
else
    netcdf.putAtt(ncid,id_skew,'standard_name','Doppler spectrum skewness at vertical polarisation');
end    
netcdf.putAtt(ncid,id_skew,'short_name','skew');
netcdf.putAtt(ncid,id_skew,'valid_range',[nanmin(data.skew(:)), nanmax(data.skew(:))]);
netcdf.putAtt(ncid,id_skew,'fill_value',str_fill_value);
netcdf.putAtt(ncid,id_skew,'comment',['Doppler spectrum skewness of measured Doppler spectrum '...
                                      'at vertical polarisation.']);

id_kurt = netcdf.defVar(ncid,'ku','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_kurt,'standard_name','RADAR.DOPPLER.SPECTRUM_KURTOSIS');
netcdf.putAtt(ncid,id_kurt,'long_name','Doppler spectrum Kurtosis');
netcdf.putAtt(ncid,id_kurt,'valid_range',[nanmin(data.kurt(:)), nanmax(data.kurt(:))]);
netcdf.putAtt(ncid,id_kurt,'fill_value','NaNf');
netcdf.putAtt(ncid,id_kurt,'comment',['Doppler spectrum Kurtosis of Doppler spectrum '...
                                    'at vertical polarisation']);                                                                 
                                  
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
    %netcdf.putAtt(ncid,id_xcorr,'valid_range',[nanmin(data.xcorr(:)), nanmax(data.xcorr(:))]);
    netcdf.putAtt(ncid,id_xcorr,'fill_value',str_fill_value);
    
    id_difphase = netcdf.defVar(ncid,'phi_dp','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_difphase,'GEOMS_name','DIFFERENTIAL.PHASE');
    netcdf.putAtt(ncid,id_difphase,'long_name','differential_phase');
    netcdf.putAtt(ncid,id_difphase,'rpg_name','difphase');    
    netcdf.putAtt(ncid,id_difphase,'standard_name','co-cross-channel differential phase');   
    netcdf.putAtt(ncid,id_difphase,'short_name','phi_dp');        
    netcdf.putAtt(ncid,id_difphase,'unite','degree');   
    %netcdf.putAtt(ncid,id_difphase,'valid_range',[nanmin(data.difphase(:)), nanmax(data.difphase(:))]); 
    netcdf.putAtt(ncid,id_difphase,'fill_value',str_fill_value); 
end


id_spec = netcdf.defVar(ncid,'spec','nc_float',[did_vel,did_range,did_time]);
netcdf.putAtt(ncid,id_spec,'GEOMS_name','RADAR.DOPPLER.SPECTRUM_VV');
if data.DualPol == 0
    netcdf.putAtt(ncid,id_spec,'standard_name','Doppler spectrum');
else
    netcdf.putAtt(ncid,id_spec,'standard_name','Doppler spectrum vertical polarization');
end        
netcdf.putAtt(ncid,id_spec,'units','mm6 mm-3');
netcdf.putAtt(ncid,id_spec,'valid_range',[nanmin(data.spec(:)), nanmax(data.spec(:))]);
netcdf.putAtt(ncid,id_spec,'fill_value','NaNf');
netcdf.putAtt(ncid,id_spec,'comment',['This is the normalized Doppler spectrum. '...
                           'The integral of the spectra minus Ze gives the noise level, ',...
                           'i.e. integral(spec)-Ze = noise. The velocity array is asymmetric, '...
                           'i.e. the absolute values of maximum and minumum velocities are not equal. ',...
                           'Since the spectrum at -v_nyquist and +v_nyquist is the same, '...
                           'the entry at +v_nyquist was cut.']);
if isfield(data, 'Ze_label') % Ze corrected, adding note
    netcdf.putAtt(ncid,id_spec,'comment2',data.Ze_label);
    netcdf.putAtt(ncid,id_spec,'corretion_dB',data.Ze_corr);
end                           

id_SLv = netcdf.defVar(ncid,'sens_limit','nc_float',[did_range,did_time]);
if data.DualPol == 0
    netcdf.putAtt(ncid,id_SLv,'standard_name','Linear sensitivity limit');
else
    netcdf.putAtt(ncid,id_SLv,'standard_name','Linear sensitivity limit for vertical polarization');
end
netcdf.putAtt(ncid,id_SLv,'units','mm6 mm-3');

id_VNoisePow_mean = netcdf.defVar(ncid,'mean_noise','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_VNoisePow_mean,'GEOMS_name','radar.doppler.spectrum_noise.power.mean_vv');
if data.DualPol == 0
    netcdf.putAtt(ncid,id_VNoisePow_mean,'standard_name','Bin Doppler spectrum mean noise power');
else
    netcdf.putAtt(ncid,id_VNoisePow_mean,'standard_name','Bin Doppler spectrum mean noise power in vertical polarization');
end    

netcdf.putAtt(ncid,id_VNoisePow_mean,'units','mm6 m-3');
netcdf.putAtt(ncid,id_VNoisePow_mean,'valid_range',[nanmin(data.VNoisePow_mean(:)), nanmax(data.VNoisePow_mean(:))]);
netcdf.putAtt(ncid,id_VNoisePow_mean,'fill_value','NaNf');
netcdf.putAtt(ncid,id_VNoisePow_mean,'comment','If cal_mom == 3, then bin noise power is only calculated if CompEna == 1.');




if data.DualPol > 0
    
    id_spec_h = netcdf.defVar(ncid,'sze_hh','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_spec_h,'standard_name','Doppler spectrum horizontal polarization');
    netcdf.putAtt(ncid,id_spec_h,'units','mm6 mm-3');
    
    id_SLh = netcdf.defVar(ncid,'sen_lim_h','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_SLh,'standard_name','Linear sensitivity limit for horizontal polarisation');
    netcdf.putAtt(ncid,id_SLh,'units','mm6 mm-3');

    if data.DualPol > 1
        id_spec_covRe = netcdf.defVar(ncid,'scov_re','nc_float',[did_vel,did_range,did_time]);
        netcdf.putAtt(ncid,id_spec_covRe,'standard_name','Real part of covariance spectrum');
        netcdf.putAtt(ncid,id_spec_covRe,'units','mm6 mm-3');
    
        id_spec_covIm = netcdf.defVar(ncid,'scov_im','nc_float',[did_vel,did_range,did_time]);
        netcdf.putAtt(ncid,id_spec_covIm,'standard_name','Imaginary part of covariance spectrum');
        netcdf.putAtt(ncid,id_spec_covIm,'units','mm6 mm-3');
    end

end % if data.DualPol > 0



if data.CompEna == 2 && data.DualPol == 2
    
    id_d_spec = netcdf.defVar(ncid,'sz_dr','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_d_spec,'standard_name','Differential spectral reflectivity');
    netcdf.putAtt(ncid,id_d_spec,'units','mm6 mm-3');
    netcdf.putAtt(ncid,id_d_spec,'valid_range',[nanmin(data.d_spec(:)), nanmax(data.d_spec(:))]);
    netcdf.putAtt(ncid,id_d_spec,'fill_value','NaNf');
    
    id_CorrCoeff = netcdf.defVar(ncid,'srho_hv','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_CorrCoeff,'standard_name','Spectral correlation coefficient');
    netcdf.putAtt(ncid,id_CorrCoeff,'valid_range',[nanmin(data.CorrCoeff(:)), nanmax(data.CorrCoeff(:))]);
    netcdf.putAtt(ncid,id_CorrCoeff,'fill_value','NaNf');
    
    id_DiffPh = netcdf.defVar(ncid,'sphi_dp','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_DiffPh,'standard_name','Spectral differential phase');
    netcdf.putAtt(ncid,id_DiffPh,'units','degree');
    netcdf.putAtt(ncid,id_DiffPh,'valid_range',[nanmin(data.DiffPh(:)), nanmax(data.DiffPh(:))]);
    netcdf.putAtt(ncid,id_DiffPh,'fill_value','NaNf');
    
end % if data.CompEna == 2 && data.DualPol > 0



if data.DualPol == 2 && data.CompEna == 2
    
    id_SLDR = netcdf.defVar(ncid,'slt_ldr','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_SLDR,'standard_name','Spectral slanted LDR');
    netcdf.putAtt(ncid,id_SLDR,'units','mm6 mm-3');
    netcdf.putAtt(ncid,id_SLDR,'valid_range',[nanmin(data.SLDR(:)), nanmax(data.SLDR(:))]);
    netcdf.putAtt(ncid,id_SLDR,'fill_value','NaNf');
    
    id_SCorrCoeff = netcdf.defVar(ncid,'sslt_ldr','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_SCorrCoeff,'standard_name','Spectral slanted corellation coefficient');
    netcdf.putAtt(ncid,id_SCorrCoeff,'valid_range',[nanmin(data.SCorrCoeff(:)), nanmax(data.SCorrCoeff(:))]);
    netcdf.putAtt(ncid,id_SCorrCoeff,'fill_value','NaNf');
    
    if data.CompEna == 2
        
        id_KDP = netcdf.defVar(ncid,'kdp','nc_float',[did_range,did_time]);
        netcdf.putAtt(ncid,id_KDP,'standard_name','Specific differential phase shift');
        netcdf.putAtt(ncid,id_KDP,'units','degree km-1');
        netcdf.putAtt(ncid,id_KDP,'valid_range',[nanmin(data.KDP(:)), nanmax(data.KDP(:))]);
        netcdf.putAtt(ncid,id_KDP,'fill_value','NaNf');
        
        id_DiffAtt = netcdf.defVar(ncid,'att_dif','nc_float',[did_range,did_time]);
        netcdf.putAtt(ncid,id_DiffAtt,'standard_name','Differential attenuation');
        netcdf.putAtt(ncid,id_DiffAtt,'units','mm6 mm-3 km-1');    
        netcdf.putAtt(ncid,id_DiffAtt,'valid_range',[nanmin(data.DiffAtt(:)), nanmax(data.DiffAtt(:))]);
        netcdf.putAtt(ncid,id_DiffAtt,'fill_value','NaNf');        
        
    end % if data.CompEna == 2
    
end % if data.DualPol == 2 && data.CompEna == 2


if data.CompEna > 0 && data.DualPol > 0
         id_HNoisePow_mean = netcdf.defVar(ncid,'h_noise_pow_mean','nc_float',[did_range,did_time]);
         netcdf.putAtt(ncid,id_HNoisePow_mean,'GEOMS_name','radar.doppler.spectrum_noise.power.mean_hh');         
         netcdf.putAtt(ncid,id_HNoisePow_mean,'standard_name','Bin Doppler spectrum mean noise power in horizontal polarization');
         netcdf.putAtt(ncid,id_HNoisePow_mean,'units','mm6 mm-3');
         netcdf.putAtt(ncid,id_HNoisePow_mean,'comment','If cal_mom == 3, then bin noise power is only calculated if CompEna == 1.');
end % data.CompEna > 0





%% ###################### initialize compression

netcdf.defVarDeflate(ncid,id_RR,true,true,9);
netcdf.defVarDeflate(ncid,id_rh,true,true,9);
netcdf.defVarDeflate(ncid,id_T_env,true,true,9);
netcdf.defVarDeflate(ncid,id_pres,true,true,9);
netcdf.defVarDeflate(ncid,id_ff,true,true,9);
netcdf.defVarDeflate(ncid,id_fff,true,true,9);
netcdf.defVarDeflate(ncid,id_dv,true,true,9);
netcdf.defVarDeflate(ncid,id_Tb,true,true,9);
if ~Flag_noLWP
    netcdf.defVarDeflate(ncid,id_lwp,true,true,9);
end

netcdf.defVarDeflate(ncid,id_vel,true,true,9);
netcdf.defVarDeflate(ncid,id_Ze,true,true,9);
netcdf.defVarDeflate(ncid,id_vm,true,true,9);
netcdf.defVarDeflate(ncid,id_sigma,true,true,9);
netcdf.defVarDeflate(ncid,id_skew,true,true,9);
netcdf.defVarDeflate(ncid,id_kurt,true,true,9);


if data.DualPol > 0
    netcdf.defVarDeflate(ncid,id_Ze_hv,true,true,9); %LP
    netcdf.defVarDeflate(ncid,id_ldr,true,true,9); %JABA    
end
if data.DualPol == 2
    netcdf.defVarDeflate(ncid,id_xcorr,true,true,9); %JABA    
    netcdf.defVarDeflate(ncid,id_difphase,true,true,9); %JABA  
end

netcdf.defVarDeflate(ncid,id_spec,true,true,9);
netcdf.defVarDeflate(ncid,id_SLv,true,true,9);
netcdf.defVarDeflate(ncid,id_VNoisePow_mean,true,true,9);


if data.DualPol > 0
    netcdf.defVarDeflate(ncid,id_spec_h,true,true,9); 
    netcdf.defVarDeflate(ncid,id_SLh,true,true,9);
    
    if data.DualPol > 1
        netcdf.defVarDeflate(ncid,id_spec_covRe,true,true,9);
        netcdf.defVarDeflate(ncid,id_spec_covIm,true,true,9);

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
netcdf.putVar(ncid,id_samples,0,data.totsamp);
netcdf.putVar(ncid,id_levels,0,data.n_levels);
netcdf.putVar(ncid,id_no_chirp_seq,0,data.no_chirp_seq);

% range dependet
netcdf.putVar(ncid,id_range,0,data.n_levels,data.range);

if exist('did_T_range', 'var')
    netcdf.putVar(ncid,id_T_levels,0,data.T_altcount,data.T_alt);
end
if exist('did_H_range', 'var')
    netcdf.putVar(ncid,id_H_levels,0,data.H_altcount,data.H_alt);
end

% chrip seq dependent variables
netcdf.putVar(ncid,id_range_offsets,0,data.no_chirp_seq,data.range_offsets);
netcdf.putVar(ncid,id_dv,0,data.no_chirp_seq,2*data.DoppMax./double(data.DoppLen));

% time dependent variables
netcdf.putVar(ncid,id_time,0,data.totsamp,data.time);
netcdf.putVar(ncid,id_RR,0,data.totsamp,data.RR);
netcdf.putVar(ncid,id_rh,0,data.totsamp,data.rh);
netcdf.putVar(ncid,id_T_env,0,data.totsamp,data.T_env);
netcdf.putVar(ncid,id_pres,0,data.totsamp,data.pres);
netcdf.putVar(ncid,id_ff,0,data.totsamp,data.ff);
netcdf.putVar(ncid,id_fff,0,data.totsamp,data.fff);
netcdf.putVar(ncid,id_Tb,0,data.totsamp,data.Tb);
if ~Flag_noLWP
    netcdf.putVar(ncid,id_lwp,0,data.totsamp,data.lwp);
end
netcdf.putVar(ncid,id_sampleTms,0,data.totsamp,data.sampleTms);


% multidimensional variables

% % % TO DO: add doppler velocity arrays
netcdf.putVar(ncid,id_vel,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.velocitymatrix,[3,2,1]));


%--- radar moments ---
netcdf.putVar(ncid,id_Ze,[0,0],[data.n_levels,data.totsamp],data.Ze');
netcdf.putVar(ncid,id_vm,[0,0],[data.n_levels,data.totsamp],data.vm');
netcdf.putVar(ncid,id_sigma,[0,0],[data.n_levels,data.totsamp],data.sigma');
netcdf.putVar(ncid,id_skew,[0,0],[data.n_levels,data.totsamp],data.skew');
netcdf.putVar(ncid,id_kurt,[0,0],[data.n_levels,data.totsamp],data.kurt');

if data.DualPol > 0
    netcdf.putVar(ncid,id_Ze_hv,[0,0],[data.n_levels,data.totsamp],data.Ze_hv'); %LP
    netcdf.putVar(ncid,id_ldr,[0,0],[data.n_levels,data.totsamp],data.LDR'); %JABA
end
if data.DualPol == 2
    disp('WARNING! skipping writing id_xcorr and difphase into files, variables not available (yet??)')
     %netcdf.putVar(ncid,id_xcorr,[0,0],[data.n_levels,data.totsamp],data.xcorr'); %JABA
     %netcdf.putVar(ncid,id_difphase,[0,0],[data.n_levels,data.totsamp],data.difphase'); %JABA
end  

netcdf.putVar(ncid,id_spec,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.spec,[3,2,1]));
netcdf.putVar(ncid,id_SLv,[0,0],[data.n_levels,data.totsamp],data.SLv');
netcdf.putVar(ncid,id_VNoisePow_mean,[0,0],[data.n_levels,data.totsamp],data.VNoisePow_mean');

if data.DualPol > 0
    netcdf.putVar(ncid,id_spec_h,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.spec_hv,[3,2,1]));
    netcdf.putVar(ncid,id_SLh,[0,0],[data.n_levels,data.totsamp],data.SLh');
    
    if data.DualPol > 1
        netcdf.putVar(ncid,id_spec_covRe,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.spec_covRe,[3,2,1]));
        netcdf.putVar(ncid,id_spec_covIm,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.spec_covIm,[3,2,1]));
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
end


netcdf.close(ncid);
