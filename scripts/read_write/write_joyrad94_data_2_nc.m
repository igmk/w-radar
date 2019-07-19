function write_joyrad94_data_2_nc(data, outfile, contactperson, processing_script)
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



%% ################ get variable ids and add attributes



%%%%%%%%%% scalar variables
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

id_AntiAlias = netcdf.defVar(ncid,'AntiAlias','nc_byte',did_scalar);
netcdf.putAtt(ncid,id_AntiAlias,'long_name','Flag for dealiasing.');
netcdf.putAtt(ncid,id_AntiAlias,'comment',...
    '0 = no dealiasing applied, 1 = dealiasing by RPG, 2 = dealiasing in process_joyrad94_data.m');

id_cal_mom = netcdf.defVar(ncid,'cal_mom','nc_byte',did_scalar);
netcdf.putAtt(ncid,id_cal_mom,'long_name','Integer indicating how moments were calculated.');
netcdf.putAtt(ncid,id_cal_mom,'comment',...
    '1 = moments were calculated from dealiased spectra. 2 = moments were calculated from raw spectra. 3 = moments were calculated by RPG software.');

id_lat = netcdf.defVar(ncid,'Lat','nc_float',did_scalar);
netcdf.putAtt(ncid,id_lat,'long_name','Latitude in degrees north [-90,90]');
netcdf.putAtt(ncid,id_lat,'units','degrees');

id_lon = netcdf.defVar(ncid,'Lon','nc_float',did_scalar);
netcdf.putAtt(ncid,id_lon,'long_name','Longitude in degrees east [-180,180]');
netcdf.putAtt(ncid,id_lon,'units','degrees');

id_MSL = netcdf.defVar(ncid,'MSL','nc_float',did_scalar);
netcdf.putAtt(ncid,id_MSL,'long_name','Height above mean sea level');
netcdf.putAtt(ncid,id_MSL,'units','m');

id_freq = netcdf.defVar(ncid,'freq','nc_float',did_scalar);
netcdf.putAtt(ncid,id_freq,'long_name','Transmission frequency');
netcdf.putAtt(ncid,id_freq,'units','GHz');


%%%%%%% range variables

id_range = netcdf.defVar(ncid,'range','nc_float',did_range);
netcdf.putAtt(ncid,id_range,'long_name','Range from antenna to the center of each range gate');
netcdf.putAtt(ncid,id_range,'units','m');

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

id_SeqAvg = netcdf.defVar(ncid,'SeqAvg','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_SeqAvg,'long_name','Number of averaged chirps in each chirp sequence');

id_SeqIntTime = netcdf.defVar(ncid,'SeqIntTime','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_SeqIntTime,'long_name','Integration time of each chirp sequence');
netcdf.putAtt(ncid,id_SeqIntTime,'units','seconds');

id_nAvg = netcdf.defVar(ncid,'nAvg','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_nAvg,'long_name','Number of spectra averaged');
netcdf.putAtt(ncid,id_nAvg,'comment','nAvg = SeqAvg/DoppLen')


id_range_offsets = netcdf.defVar(ncid,'range_offsets','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_range_offsets,'long_name','Chirp sequence start index array in range array');
netcdf.putAtt(ncid,id_range_offsets,'comment',...
    'The command range(range_offsets) will give you the range where a new chirp sequence starts. range_offsets counts from 1 to n_levels.');

id_dr = netcdf.defVar(ncid,'dr','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_dr,'long_name','Range resolution for chirp sequences');
netcdf.putAtt(ncid,id_dr,'units','m');

id_DoppLen = netcdf.defVar(ncid,'DoppLen','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_DoppLen,'long_name','Number of samples in Dopppler spectra of each chirp sequence. Needed to calculate the Doppler resolution: DoppRes = 2*DoppMax/DoppLen');

id_DoppMax = netcdf.defVar(ncid,'DoppMax','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_DoppMax,'long_name','Max. unambigious Doppler velocity for each chirp sequence. Needed to calculate the Doppler resolution: DoppRes = 2*DoppMax/DoppLen');
netcdf.putAtt(ncid,id_DoppMax,'units','m/s');





%%%%%%%% time dependend variables

id_time = netcdf.defVar(ncid,'time','nc_uint',did_time);
netcdf.putAtt(ncid,id_time,'long_name','Time in sec since 2001.01.01. 00:00:00');
netcdf.putAtt(ncid,id_time,'units','seconds UTC');
netcdf.putAtt(ncid,id_time,'comment','To get the correct time the variable sampleTms must be added: time = time + sampleTms.');

id_sampleTms = netcdf.defVar(ncid,'sampleTms','nc_int',did_time);
netcdf.putAtt(ncid,id_sampleTms,'long_name','Milliseconds of sample');
netcdf.putAtt(ncid,id_sampleTms,'units','ms');
netcdf.putAtt(ncid,id_sampleTms,'comment','To get the correct time the variable sampleTms must be added: time = time + sampleTms.');

id_RR = netcdf.defVar(ncid,'RR','nc_float',did_time);
netcdf.putAtt(ncid,id_RR,'long_name','Rain rate of meteo-station');
netcdf.putAtt(ncid,id_RR,'units','mm/h');

id_rh = netcdf.defVar(ncid,'rh','nc_float',did_time);
netcdf.putAtt(ncid,id_rh,'long_name','Relative humidity of meteo-station');
netcdf.putAtt(ncid,id_rh,'units','%');

id_T_env = netcdf.defVar(ncid,'T_env','nc_float',did_time);
netcdf.putAtt(ncid,id_T_env,'long_name','Environmental temperature of meteo-station');
netcdf.putAtt(ncid,id_T_env,'units','K');

id_pres = netcdf.defVar(ncid,'pres','nc_float',did_time);
netcdf.putAtt(ncid,id_pres,'long_name','Environmental pressure of meteo-station');
netcdf.putAtt(ncid,id_pres,'units','hPa');

id_ff = netcdf.defVar(ncid,'ff','nc_float',did_time);
netcdf.putAtt(ncid,id_ff,'long_name','Wind speed of meteo-station');
netcdf.putAtt(ncid,id_ff,'units','km/h');

id_fff = netcdf.defVar(ncid,'fff','nc_float',did_time);
netcdf.putAtt(ncid,id_fff,'long_name','Wind direction of meteo-station');
netcdf.putAtt(ncid,id_fff,'units','degrees');

id_vol = netcdf.defVar(ncid,'vol','nc_float',did_time);
netcdf.putAtt(ncid,id_vol,'long_name','Voltage direct detection channel');
netcdf.putAtt(ncid,id_vol,'units','V');

id_Tb = netcdf.defVar(ncid,'Tb','nc_float',did_time);
netcdf.putAtt(ncid,id_Tb,'long_name','brightness temperature direct detection channel');
netcdf.putAtt(ncid,id_Tb,'units','K');

id_lwp = netcdf.defVar(ncid,'lwp','nc_float',did_time);
netcdf.putAtt(ncid,id_lwp,'long_name','Liquid water path calculated by RPG software');
netcdf.putAtt(ncid,id_lwp,'units','g/m^2');

id_powIF = netcdf.defVar(ncid,'PowIF','nc_float',did_time);
netcdf.putAtt(ncid,id_powIF,'long_name','IF power at ADC');
netcdf.putAtt(ncid,id_powIF,'units','microWatts');

id_ele = netcdf.defVar(ncid,'ele','nc_float',did_time);
netcdf.putAtt(ncid,id_ele,'long_name','Elevation angle');
netcdf.putAtt(ncid,id_ele,'units','degrees');

id_az = netcdf.defVar(ncid,'az','nc_float',did_time);
netcdf.putAtt(ncid,id_az,'long_name','Azimuth angle');
netcdf.putAtt(ncid,id_az,'units','degrees');

id_status = netcdf.defVar(ncid,'status','nc_float',did_time);
netcdf.putAtt(ncid,id_status,'long_name','status flag: 0/1 = heater on/off; 0/10 = blower on/off');

id_TransPow = netcdf.defVar(ncid,'TransPow','nc_float',did_time);
netcdf.putAtt(ncid,id_TransPow,'long_name','Transmitted power');
netcdf.putAtt(ncid,id_TransPow,'units','W');

id_T_trans = netcdf.defVar(ncid,'T_trans','nc_float',did_time);
netcdf.putAtt(ncid,id_T_trans,'long_name','Transmitter temperature');
netcdf.putAtt(ncid,id_T_trans,'units','K');

id_T_rec = netcdf.defVar(ncid,'T_rec','nc_float',did_time);
netcdf.putAtt(ncid,id_T_rec,'long_name','Receiver temperature');
netcdf.putAtt(ncid,id_T_rec,'units','K');

id_T_pc = netcdf.defVar(ncid,'T_pc','nc_float',did_time);
netcdf.putAtt(ncid,id_T_pc,'long_name','PC temperature');
netcdf.putAtt(ncid,id_T_pc,'units','K');

id_QF = netcdf.defVar(ncid,'QF','nc_byte',did_time);
netcdf.putAtt(ncid,id_QF,'long_name','Quality flag given by radar');
netcdf.putAtt(ncid,id_QF,'comment', ...
    ['To get the bit entries, one has to convert the integer into a 4 bit binary. '...
    'bit4 = ADC saturation, bit3 = spectral width too high, bit2 = no transmitter power leveling.' ...
    'Note that in the above convention holds: bit1 = 2^3, bit2 = 2^2, bit3 = 2^1, bit4 = 2^0'])





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

id_Aliasmask = netcdf.defVar(ncid,'AliasMask','nc_byte',[did_range,did_time]);
netcdf.putAtt(ncid,id_Aliasmask,'long_name','Mask array indicating in which bin dealiasing was applied. If AnitAlias = 1, then dealiasing was applied by RPG software: 0 = not applied; 1 = applied; If AntiAlias = 2, then dealiasing was applied in post-processing: 0 = no aliasing detected, 1 = aliasing detected; if any bin equals 1 (while AntiAlias = 2) then the full column was dealiased.');

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

id_Ze = netcdf.defVar(ncid,'Ze','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_Ze,'long_name','Equivalent radar reflectivity factor Ze');
netcdf.putAtt(ncid,id_Ze,'units','mm^6/m^3');

if isfield(data, 'Ze_label') % Ze corrected, adding note
    netcdf.putAtt(ncid,id_Ze,'comment',data.Ze_label);
    netcdf.putAtt(ncid,id_Ze,'corretion_dB',data.Ze_corr);
    
end

id_vm = netcdf.defVar(ncid,'vm','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_vm,'long_name','Mean Doppler velocity');
netcdf.putAtt(ncid,id_vm,'units','m/s');
netcdf.putAtt(ncid,id_vm,'comment','negative values indicate falling particles towards the radar')

id_sigma = netcdf.defVar(ncid,'sigma','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_sigma,'long_name','Spectral width of Doppler velocity spectrum');
netcdf.putAtt(ncid,id_sigma,'units','m/s');

id_skew = netcdf.defVar(ncid,'skew','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_skew,'long_name','Skewness');

id_kurt = netcdf.defVar(ncid,'kurt','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_kurt,'long_name','Kurtosis');

id_PNv = netcdf.defVar(ncid,'PNv','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_PNv,'long_name','Total IF power in vertical polarization measured at ADC input');
netcdf.putAtt(ncid,id_PNv,'units','W');

id_SLv = netcdf.defVar(ncid,'SLv','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_SLv,'long_name','Linear sensitivity limit for vertical polarisation');
netcdf.putAtt(ncid,id_SLv,'units','mm^6/mm^3');

id_spec = netcdf.defVar(ncid,'spec','nc_float',[did_vel,did_range,did_time]);
netcdf.putAtt(ncid,id_spec,'long_name','Doppler spectrum vertical polarization');
netcdf.putAtt(ncid,id_spec,'units','mm^6/mm^3');
netcdf.putAtt(ncid,id_spec,'comment',['This is the normalized Doppler spectrum. The integral of the spectra minus Ze gives the noise level, ',...
    'i.e. integral(spec)-Ze = noise. The velocity array is asymmetric, i.e. the absolute values of maximum and minumum velocities are not equal. ',...
    'Since the spectrum at -v_nyquist and +v_nyquist is the same, the entry at +v_nyquist was cut.']);

id_VNoisePow_mean = netcdf.defVar(ncid,'VNoisePow_mean','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_VNoisePow_mean,'long_name','Bin Doppler spectrum mean noise power in vertical polarization');
netcdf.putAtt(ncid,id_VNoisePow_mean,'units','mm^6/m^3');
netcdf.putAtt(ncid,id_VNoisePow_mean,'comment','If cal_mom == 3, then bin noise power is only calculated if CompEna == 1.');

id_VNoisePow_peak = netcdf.defVar(ncid,'VNoisePow_peak','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_VNoisePow_peak,'long_name','Bin Doppler spectrum peak noise power in vertical polarization');
netcdf.putAtt(ncid,id_VNoisePow_peak,'units','mm^6/m^3');
netcdf.putAtt(ncid,id_VNoisePow_peak,'comment','If cal_mom == 3, then bin noise power is only calculated if CompEna == 1.');

id_QualFlag = netcdf.defVar(ncid,'QualityFlag','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_QualFlag,'long_name','Quality flag, added in the additional data processing to alert for known issues');
netcdf.putAtt(ncid,id_QualFlag,'comment', ...
    ['This variable contains information on anything that might impact the quality ', ...
    'of the data at each pixel. Must be converted into three bit binary string. ', ...
	'If 0, i.e. dec2bin(QualityFlag,3) = 000, none of the included issues were ', ...
    'found. The definitions of each bit are given in the definition attribute.']);
netcdf.putAtt(ncid,id_QualFlag,'definition', ...
    ['If 2^0 bit is 1: this range gate is known to have aritifical spikes occurring', ...
     'If 2^1 bit is 1: aircraft or other known flying non-meteorological object', ...
     'If 2^2 bit is 1: wet-radome (was a problem for mirac-a for a time period with coating missing)' ...
     ]);



if data.DualPol > 0
    
    id_PNh = netcdf.defVar(ncid,'PNh','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_PNh,'long_name','Total IF power in horizontal polarization measured at ADC input');
    netcdf.putAtt(ncid,id_PNh,'units','W');
    
    id_SLh = netcdf.defVar(ncid,'SLh','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_SLh,'long_name','Linear sensitivity limit for vertical polarisation');
    netcdf.putAtt(ncid,id_SLh,'units','mm^6/mm^3');
    
    id_spec_h = netcdf.defVar(ncid,'spec_h','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_spec_h,'long_name','Doppler spectrum horizontal polarization');
    netcdf.putAtt(ncid,id_spec_h,'units','mm^6/mm^3');
    
    id_spec_covRe = netcdf.defVar(ncid,'spec_covRe','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_spec_covRe,'long_name','Real part of covariance spectrum');
    netcdf.putAtt(ncid,id_spec_covRe,'units','mm^6/mm^3)');
    
    id_spec_covIm = netcdf.defVar(ncid,'spec_covIm','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_spec_covIm,'long_name','Real part of covariance spectrum');
    netcdf.putAtt(ncid,id_spec_covIm,'units','mm^6/mm^3)');

    %Included by Bravo-Aranda, J.A. JABA
    id_ldr = netcdf.defVar(ncid,'ldr','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_ldr,'long_name','Linear depolarization ratio');
    netcdf.putAtt(ncid,id_ldr,'units','');    
    
    if data.DualPol == 2
        id_xcorr = netcdf.defVar(ncid,'xcorr','nc_float',[did_range,did_time]);
        netcdf.putAtt(ncid,id_xcorr,'long_name','co-cross-channel correlation coefficient');
        netcdf.putAtt(ncid,id_xcorr,'units','');    

        id_difphase = netcdf.defVar(ncid,'difphase','nc_float',[did_range,did_time]);
        netcdf.putAtt(ncid,id_difphase,'long_name','co-cross-channel differential phase');
        netcdf.putAtt(ncid,id_difphase,'units','');        
    end
end % if data.DualPol > 0


if data.CompEna == 2 && data.DualPol == 2
    
    id_d_spec = netcdf.defVar(ncid,'d_spec','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_d_spec,'long_name','Differntial spectral reflectivity');
    netcdf.putAtt(ncid,id_d_spec,'units','dB');
    
    id_CorrCoeff = netcdf.defVar(ncid,'CorrCoeff','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_CorrCoeff,'long_name','Spectral correlation coefficient');
    
    id_DiffPh = netcdf.defVar(ncid,'DiffPh','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_DiffPh,'long_name','Spectral differential phase');
    netcdf.putAtt(ncid,id_DiffPh,'units','dB');
    
end % if data.CompEna == 2 && data.DualPol > 0


if data.DualPol == 2 && data.CompEna == 2
    
    id_SLDR = netcdf.defVar(ncid,'SLDR','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_SLDR,'long_name','Spectral slanted LDR');
    netcdf.putAtt(ncid,id_SLDR,'units','dB');
    
    id_SCorrCoeff = netcdf.defVar(ncid,'SCOrrCoeff','nc_float',[did_vel,did_range,did_time]);
    netcdf.putAtt(ncid,id_SCorrCoeff,'long_name','Spectral slanted corellation coefficient');
    
    if data.CompEna == 2
        
        id_KDP = netcdf.defVar(ncid,'KDP','nc_float',[did_range,did_time]);
        netcdf.putAtt(ncid,id_KDP,'long_name','Specific differential phase shift');
        netcdf.putAtt(ncid,id_KDP,'units','rad/km');
        
        id_DiffAtt = netcdf.defVar(ncid,'DiffAtt','nc_float',[did_range,did_time]);
        netcdf.putAtt(ncid,id_DiffAtt,'long_name','Differential attenuation');
        netcdf.putAtt(ncid,id_DiffAtt,'units','dB/km');        
        
    end % if data.CompEna == 2
    
end % if data.DualPol == 2 && data.CompEna == 2

if data.CompEna > 0 && data.DualPol > 0
         id_HNoisePow_mean = netcdf.defVar(ncid,'HNoisePow_mean','nc_float',[did_range,did_time]);
         netcdf.putAtt(ncid,id_HNoisePow_mean,'long_name','Bin Doppler spectrum mean noise power in horizontal polarization');
         netcdf.putAtt(ncid,id_HNoisePow_mean,'units','mm^6/m^3');
         netcdf.putAtt(ncid,id_HNoisePow_mean,'comment','If cal_mom == 3, then bin noise power is only calculated if CompEna == 1.');
         
         id_HNoisePow_peak = netcdf.defVar(ncid,'HNoisePow_peak','nc_float',[did_range,did_time]);
         netcdf.putAtt(ncid,id_HNoisePow_peak,'long_name','Bin Doppler spectrum peak noise power in horizontal polarization');
         netcdf.putAtt(ncid,id_HNoisePow_peak,'units','mm^6/m^3');
         netcdf.putAtt(ncid,id_HNoisePow_peak,'comment','If cal_mom == 3, then bin noise power is only calculated if CompEna == 1.');
end % data.CompEna > 0


%% ######################## add global attributes
glob = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,glob,'FillValue','NaN');
netcdf.putAtt(ncid,glob,'program_name',data.progname);
if data.modelno == 0
    model = '94 GHz single pol.';
else
    model = '94 GHz dual pol.';
end
netcdf.putAtt(ncid,glob,'model_type',model);
netcdf.putAtt(ncid,glob,'contact',contactperson);
netcdf.putAtt(ncid,glob,'processing script',processing_script);


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
    netcdf.putVar(ncid,id_spec_h,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.spec_h,[3,2,1]));
    netcdf.putVar(ncid,id_spec_covRe,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.spec_covRe,[3,2,1]));
    netcdf.putVar(ncid,id_spec_covIm,[0,0,0],[max(data.DoppLen),data.n_levels,data.totsamp],permute(data.spec_covIm,[3,2,1]));
    netcdf.putVar(ncid,id_ldr,[0,0],[data.n_levels,data.totsamp],data.LDR'); %JABA
    if data.DualPol == 2
        netcdf.putVar(ncid,id_xcorr,[0,0],[data.n_levels,data.totsamp],data.xcorr'); %JABA
        netcdf.putVar(ncid,id_difphase,[0,0],[data.n_levels,data.totsamp],data.difphase'); %JABA
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
