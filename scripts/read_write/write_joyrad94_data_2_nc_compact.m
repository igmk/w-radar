function write_joyrad94_data_2_nc_compact(data,outfile, config)

% this function writes joyrad94 data into netcdf4

%% ################## Create a netCDF file.

ncid = netcdf.create(outfile,'NETCDF4'); 


%% ################# Define dimensions

%did_time = netcdf.defDim(ncid,'time',data.totsamp);
%Changed by J.A. Bravo-Aranda
did_time = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
did_range = netcdf.defDim(ncid,'range',data.n_levels);
did_no_seq = netcdf.defDim(ncid,'chirp_sequences',data.no_chirp_seq);
did_scalar = netcdf.defDim(ncid,'scalar',1);


%% ################ get variable ids and add attributes

%%%%%%%%%% scalar variables

id_AntiAlias = netcdf.defVar(ncid,'AntiAlias','nc_byte',did_scalar);
netcdf.putAtt(ncid,id_AntiAlias,'long_name','Flag for dealiasing.');
netcdf.putAtt(ncid,id_AntiAlias,'comment',...
    '0 = no dealiasing applied, 1 = dealiasing by RPG, 2 = dealiasing in process_joyrad94_data.m');

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



%%%%%%%% chirp_seq_dependent variables

id_SeqAvg = netcdf.defVar(ncid,'SeqAvg','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_SeqAvg,'long_name','Number of averaged chirps in each chirp sequence');

id_SeqIntTime = netcdf.defVar(ncid,'SeqIntTime','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_SeqIntTime,'long_name','Integration time of each chirp sequence');
netcdf.putAtt(ncid,id_SeqIntTime,'units','seconds');

id_DoppLen = netcdf.defVar(ncid,'DoppLen','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_DoppLen,'long_name','Number of samples in Dopppler spectra of each chirp sequence. Needed to calculate the Doppler resolution: DoppRes = 2*DoppMax/DoppLen');

id_DoppMax = netcdf.defVar(ncid,'DoppMax','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_DoppMax,'long_name','Max. unambigious Doppler velocity for each chirp sequence. Needed to calculate the Doppler resolution: DoppRes = 2*DoppMax/DoppLen');
netcdf.putAtt(ncid,id_DoppMax,'units','m/s');

id_nAvg = netcdf.defVar(ncid,'nAvg','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_nAvg,'long_name','Number of spectra averaged');
netcdf.putAtt(ncid,id_nAvg,'comment','nAvg = SeqAvg/DoppLen')


id_range_offsets = netcdf.defVar(ncid,'range_offsets','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_range_offsets,'long_name','Chirp sequence start index array in range array');
netcdf.putAtt(ncid,id_range_offsets,'comment',...
    'The command range(range_offsets) will give you the range where a new chirp sequence starts. range_offsets counts from 1 to n_levels.');




%%%%%%%% time dependend variables

id_time = netcdf.defVar(ncid,'time','nc_uint',did_time);
netcdf.putAtt(ncid,id_time,'long_name','Time in sec since 2001.01.01. 00:00:00');
netcdf.putAtt(ncid,id_time,'units','seconds UTC');
netcdf.putAtt(ncid,id_time,'comment','To get the correct time the variable sampleTms must be added: time = time + sampleTms.');
if isfield( data, totsampchangelabel )
    netcdf.putAtt(ncid,id_time, 'Flag', 'Dublicate time stamps found in lv0-file, the first occurrence of the dublicate time is removed')
end


id_sampleTms = netcdf.defVar(ncid,'sampleTms','nc_int',did_time);
netcdf.putAtt(ncid,id_sampleTms,'long_name','Milliseconds of sample');
netcdf.putAtt(ncid,id_sampleTms,'units','mu s');
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

id_Tb = netcdf.defVar(ncid,'Tb','nc_float',did_time);
netcdf.putAtt(ncid,id_Tb,'long_name','brightness temperature direct detection channel');
netcdf.putAtt(ncid,id_Tb,'units','K');

id_lwp = netcdf.defVar(ncid,'lwp','nc_float',did_time);
netcdf.putAtt(ncid,id_lwp,'long_name','Liquid water path calculated by RPG software');
netcdf.putAtt(ncid,id_lwp,'units','g/m^2');


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
     'If 2^2 bit is 1: wet-radome (was a problem for mirac-a for a time period when coating missing)' ...
     ]);

id_Aliasmask = netcdf.defVar(ncid,'AliasMask','nc_byte',[did_range,did_time]);
netcdf.putAtt(ncid,id_Aliasmask,'long_name','Mask array indicating in which bin dealiasing was applied. If AnitAlias = 1, then dealiasing was applied by RPG software: 0 = not applied; 1 = applied; If AntiAlias = 2, then dealiasing was applied in post-processing: 0 = no aliasing detected, 1 = aliasing detected; if any bin equals 1 (while AntiAlias = 2) then the full column was dealiased.');

%Included by Bravo-Aranda, J.A. JABA
id_ldr = netcdf.defVar(ncid,'ldr','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_ldr,'long_name','Linear depolarization ratio');

if data.DualPol == 2
    id_xcorr = netcdf.defVar(ncid,'xcorr','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_xcorr,'long_name','co-cross-channel correlation coefficient');

    id_difphase = netcdf.defVar(ncid,'difphase','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_difphase,'long_name','co-cross-channel differential phase');
    
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
netcdf.putAtt(ncid,glob,'model_type',model);
netcdf.putAtt(ncid,glob,'contact',config.contactperson);
netcdf.putAtt(ncid,glob,'processing script',config.processing_script);


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
netcdf.putVar(ncid,id_freq,0,data.freq);
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
netcdf.putVar(ncid,id_ff,0,data.totsamp,data.ff);
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
    netcdf.putVar(ncid,id_xcorr,[0,0],[data.n_levels,data.totsamp],data.xcorr'); %JABA
    netcdf.putVar(ncid,id_difphase,[0,0],[data.n_levels,data.totsamp],data.difphase'); %JABA
    
end

netcdf.close(ncid);

end
