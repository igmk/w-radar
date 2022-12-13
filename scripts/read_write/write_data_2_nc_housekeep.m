function write_data_2_nc_housekeep(data, outfile, config)
% function to write house keeping data and any other variable and parameter
% not included in moments or spectra file (incl. 89 GHz TB and weather
% station data)
% RG 1.6.2022

%################## Create a netCDF file.

ncid = netcdf.create(outfile,'NETCDF4'); 

%% ################# Define dimensions
% did_time   = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
did_time   = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
did_range  = netcdf.defDim(ncid,'range',data.n_levels);
did_no_seq = netcdf.defDim(ncid,'chirp_sequence',data.no_chirp_seq);
did_str = netcdf.defDim(ncid,'string',50); % assuming 50 charachters would be enough

Temp_ret_flag = false;
if isfield(data, 'T_altcount')
    if ne(data.T_altcount,0)
        did_T_range = netcdf.defDim(ncid,'T_levels',data.T_altcount);
        Temp_ret_flag = true;
    end    
end

Hum_ret_flag = false;
if isfield(data, 'H_altcount')
    if ne(data.H_altcount,0)
        did_H_range = netcdf.defDim(ncid,'H_levels',data.H_altcount);
        Hum_ret_flag = true;
    end
end


%% ######################## add global attributes

write_globatt_ac3(ncid, config, data, 3)

if data.cal_mom == 3 % if moments not calculated by our script, add note
    netcdf.putAtt(ncid,glob,'note_dataprocessing', 'Moments are taken from RPG processing sofware (lv1)');
end


%% ################ get variable ids and add attributes
% collection of functions defining variables with attributes
defh = outvarmeta;

%%%%%%%%%% coordinate variables
id_time = defh.time(ncid, did_time, 'nc_int', isfield(data, 'totsampchangelabel' ));
netcdf.putAtt(ncid,id_time,'rpg_name','Time');

id_range = defh.range(ncid, did_range);
netcdf.putAtt(ncid,id_range,'rpg_name','RAlts');

id_no_seq = defh.chirp_sequence(ncid, did_no_seq);

if Temp_ret_flag
    id_coord_Tlevels = netcdf.defVar(ncid,'T_levels','nc_float',did_T_range);
    netcdf.putAtt(ncid,id_coord_Tlevels,'long_name','height levels for retrieved temperature profile');
    netcdf.putAtt(ncid,id_coord_Tlevels,'units','m');
    netcdf.putAtt(ncid,id_coord_Tlevels,'positive','up'); % recommended by CF convention
    netcdf.putAtt(ncid,id_coord_Tlevels,'rpg_name','TAlts');
end

if Hum_ret_flag
    id_coord_Hlevels = netcdf.defVar(ncid,'H_levels','nc_float',did_H_range);
    netcdf.putAtt(ncid,id_coord_Hlevels,'long_name','height levels for retrieved humidity profiles');
    netcdf.putAtt(ncid,id_coord_Hlevels,'units','m');
    netcdf.putAtt(ncid,id_coord_Hlevels,'positive','up'); % recommended by CF convention
    netcdf.putAtt(ncid,id_coord_Hlevels,'rpg_name','HAlts');
end

id_sampleTms = defh.sampleTms(ncid, did_time);
netcdf.putAtt(ncid,id_sampleTms,'rpg_name','MSec');

id_height = defh.height(ncid, did_range);
netcdf.defVarFill(ncid,id_height,false,NaN('single'))


%%%%%%%%%% scalar variables

id_filecode = netcdf.defVar(ncid,'filecode','nc_int',[]);
netcdf.putAtt(ncid,id_filecode,'long_name','filecode/version number of lv0 file');
netcdf.defVarFill(ncid,id_filecode,false,int32(-999))
netcdf.putAtt(ncid,id_filecode,'rpg_name','File Code'); 


id_progno = netcdf.defVar(ncid,'chirp_program_number','nc_int',[]);
netcdf.putAtt(ncid,id_progno,'long_name','chirp program number');
netcdf.defVarFill(ncid,id_progno,false,int32(-999))
netcdf.putAtt(ncid,id_progno,'comment','used to identify specific chirp tables');
netcdf.putAtt(ncid,id_progno,'rpg_name','ProgNo');

if isfield(data, 'custname')
    id_custname = netcdf.defVar(ncid,'customer_name','nc_char',did_str);
    netcdf.putAtt(ncid,id_custname,'comment','customer name as given in RPG software');
    netcdf.putAtt(ncid,id_custname,'rpg_name','CustName');
end

% commented out, because does not provide any info - the variables range,
% time and chirp_sequence alrady tells how many bins in the data. RG 20.9.2022
% 
% id_samples = netcdf.defVar(ncid,'n_profiles','nc_int',[]);
% netcdf.putAtt(ncid,id_samples,'long_name','Number of profiles');
% netcdf.putAtt(ncid,id_samples,'rpg_name','TotSamp'); 
% 
% id_levels = netcdf.defVar(ncid,'n_levels','nc_int',[]);
% netcdf.putAtt(ncid,id_levels,'long_name','Number of range gates');
% netcdf.putAtt(ncid,id_levels,'rpg_name','RAltN'); 
% 
% id_no_chirp_seq = netcdf.defVar(ncid,'no_chirp_seq','nc_int',[]);
% netcdf.putAtt(ncid,id_no_chirp_seq,'long_name','Number of chirp sequences');
% netcdf.putAtt(ncid,id_no_chirp_seq,'comment',['The radar can be programmed '...
%                                               'to run different resolution ' ...
%                                               'modes for different layers ' ...
%                                               'spanning several range gates']);


id_AntiAlias = netcdf.defVar(ncid,'alias_method_flag','nc_byte',[]);
netcdf.putAtt(ncid,id_AntiAlias,'long_name','Flag for applied dealiasing method');
netcdf.putAtt(ncid,id_AntiAlias,'comment',['The flag index shows: '...
                                           '0 = no dealiasing applied, '... 
                                           '1 = dealiasing by RPG, '...
                                           '2 = dealiasing by processing_script']);


id_pol = netcdf.defVar(ncid,'DualPol','nc_byte',[]);
netcdf.putAtt(ncid,id_pol,'long_name','Polarisation status')
netcdf.putAtt(ncid,id_pol,'comment','DualPol = 2: Full polarimetry. DualPol = 1: Polarimetric receiver but only vertical polarization transmitted. DualPol = 0: Only vertical polarization available.');
netcdf.putAtt(ncid,id_pol,'rpg_name','DualPol'); 

id_T_levels = netcdf.defVar(ncid,'n_T_levels','nc_int',[]);
netcdf.putAtt(ncid,id_T_levels,'long_name','Number of range levels of retrieved temp. profile');
netcdf.defVarFill(ncid,id_T_levels,false,int32(-999))
netcdf.putAtt(ncid,id_T_levels,'rpg_name','TAltN'); 

id_H_levels = netcdf.defVar(ncid,'n_H_levels','nc_int',[]);
netcdf.putAtt(ncid,id_H_levels,'long_name','Number of range levels of retrieved hum. profile');
netcdf.defVarFill(ncid,id_H_levels,false,int32(-999))
netcdf.putAtt(ncid,id_H_levels,'rpg_name','HAltN'); 
                                       
id_compress = netcdf.defVar(ncid,'compression_flag','nc_byte',[]);
netcdf.putAtt(ncid,id_compress,'long_name','flag for spectra compression in RPG software');
netcdf.putAtt(ncid,id_compress,'compression','If compression_flag = 1, then compression was enabled in RPG software and noise has been removed from spectra, otherwise 0.');
if ~isfield(data, 'CompEna')
    data.CompEna = 0;
end

% not a good way to include this information
% id_cal_mom = netcdf.defVar(ncid,'cal_mom','nc_byte',[]);
% netcdf.putAtt(ncid,id_cal_mom,'long_name','Integer indicating how moments were calculated.');
% netcdf.putAtt(ncid,id_cal_mom,'comment',['Set in config file. ' ...
%                                          '2 = moments were calculated from ' ...
%                                          'spectra. 3 = moments are taken ' ...
%                                          'from RPG processing sofware (lv1).']);

id_swv = netcdf.defVar(ncid,'radar_software','nc_char',did_str);
netcdf.putAtt(ncid,id_swv,'long_name','radar software version');


id_AntSep = netcdf.defVar(ncid,'antenna_separation','nc_float',[]);
netcdf.putAtt(ncid,id_AntSep,'long_name','Separation of antenna axis');
netcdf.putAtt(ncid,id_AntSep,'units','m');
netcdf.defVarFill(ncid,id_AntSep,false,NaN('single'))
netcdf.putAtt(ncid,id_AntSep,'rpg_name','AntSep'); 

id_AntDia = netcdf.defVar(ncid,'antenna_diameter','nc_float',[]);
netcdf.putAtt(ncid,id_AntDia,'long_name','Antenna diameter');
netcdf.putAtt(ncid,id_AntDia,'units','m');
netcdf.defVarFill(ncid,id_AntDia,false,NaN('single'))
netcdf.putAtt(ncid,id_AntDia,'rpg_name','AntDia'); 
 
if strcmp(data.radarsw, '1.00') % the "radar constant" variable given in data for first software version is defined differently (see manual for details)
    id_C = netcdf.defVar(ncid,'radar_constant','nc_float',did_time);
    netcdf.putAtt(ncid,id_C,'long_name','Radar constant');
    netcdf.defVarFill(ncid,id_C,false,NaN('single'))
    netcdf.putAtt(ncid,id_C,'comment','note that radar constant given in different format for software version 1, see manual for details');
    netcdf.putAtt(ncid,id_C,'rpg_name','RadC');
else
    id_C = netcdf.defVar(ncid,'radar_constant','nc_float',[]);
    netcdf.putAtt(ncid,id_C,'long_name','Radar constant');
    netcdf.defVarFill(ncid,id_C,false,NaN('single'))
    netcdf.putAtt(ncid,id_C,'comment','see manual for details');
    netcdf.putAtt(ncid,id_C,'rpg_name','Cr'); 
end

if isfield(data, 'SampRate')
    id_smprt = netcdf.defVar(ncid,'ADC_sampling_rate','nc_float',[]);
    netcdf.putAtt(ncid,id_smprt,'comment','ADC sampling rate');
    netcdf.putAtt(ncid,id_smprt,'units','Hz');
    netcdf.defVarFill(ncid,id_smprt,false,NaN('single'))
    netcdf.putAtt(ncid,id_smprt,'rpg_name','SampRate'); % data.SampRate
end

if isfield(data, 'MaxRange')
    id_maxr = netcdf.defVar(ncid,'max_range','nc_float',[]);
    netcdf.putAtt(ncid,id_maxr,'long_name','maximum unambiguous range');
    netcdf.putAtt(ncid,id_maxr,'units','m');
    netcdf.defVarFill(ncid,id_maxr,false,NaN('single'))
    netcdf.putAtt(ncid,id_maxr,'rpg_name','MaxRange'); % data.MaxRange
end

if isfield(data, 'SupPowLev')
    id_powesup = netcdf.defVar(ncid,'power_level_suppression','nc_float',[]);
    netcdf.putAtt(ncid,id_powesup,'long_name','Power level suppression');
    netcdf.defVarFill(ncid,id_powesup,false,NaN('single'))
    netcdf.putAtt(ncid,id_powesup,'comment','indicates if power levelling has been used (0=yes, 1=no)');
    netcdf.putAtt(ncid,id_powesup,'rpg_name','SupPowLev'); % data.SupPowLev
end

if isfield(data, 'SpkFilEna')
    id_sfil = netcdf.defVar(ncid,'spike_filter','nc_float',[]);
    netcdf.putAtt(ncid,id_sfil,'long_name','spike filter');
    netcdf.defVarFill(ncid,id_sfil,false,NaN('single'))
    netcdf.putAtt(ncid,id_sfil,'comment','indicates if spike/plankton filter has been used in RPG software(0=yes, 1=no)');
    netcdf.putAtt(ncid,id_sfil,'rpg_name','SpkFilEna'); % data.SpkFilEna
end

if isfield(data, 'FFTWindow')
    id_fftwin = netcdf.defVar(ncid,'FFT_window','nc_char',did_str);
    netcdf.putAtt(ncid,id_fftwin,'long_name','FFT window');
    netcdf.putAtt(ncid,id_fftwin,'comment','FFT window used in signal processing');
    netcdf.putAtt(ncid,id_fftwin,'rpg_name','FFTWindow'); % data.FFTWindow
end

if isfield(data, 'FFTInputRng')
    id_finrn = netcdf.defVar(ncid,'ADC_input_voltage_range','nc_float',[]);
    netcdf.putAtt(ncid,id_finrn,'long_name','ADC input voltage range (+/-)');
    netcdf.putAtt(ncid,id_finrn,'units','mV');
    netcdf.defVarFill(ncid,id_finrn,false,NaN('single'))
    netcdf.putAtt(ncid,id_finrn,'rpg_name','FFTInputRng'); % data.FFTInputRng
end

if isfield(data, 'NoiseFilt')
    id_noisfil = netcdf.defVar(ncid,'noise_filter','nc_float',[]);
    netcdf.putAtt(ncid,id_noisfil,'long_name','noise filter factor');
    netcdf.defVarFill(ncid,id_noisfil,false,NaN('single'))
    netcdf.putAtt(ncid,id_noisfil,'comment','Determines the threshold used for data logging, which is the STD of noise multiplied with the noise filter factor');
    netcdf.putAtt(ncid,id_noisfil,'rpg_name','NoiseFilt'); % data.NoiseFilt
end

id_SampDur = netcdf.defVar(ncid,'sample_duration','nc_float',[]);
netcdf.putAtt(ncid,id_SampDur,'long_name','full sample duration');
netcdf.putAtt(ncid,id_SampDur,'units','sec');
netcdf.defVarFill(ncid,id_SampDur,false,NaN('single'))
netcdf.putAtt(ncid,id_SampDur,'comment','This is the sum of the chirp sequences integration times, and not the integration time (see chirp_seq_int_time)')


%%%% for polarimetric radar
if data.DualPol > 0

    id_phcor = netcdf.defVar(ncid,'phase_correction','nc_float',[]);
    netcdf.putAtt(ncid,id_phcor,'comment','flag indicating, if phase correction (dual. pol. radars) has been used (0=yes, 1=no)');
    netcdf.putAtt(ncid,id_phcor,'rpg_name','PhaseCorr'); % data.PhaseCorr

    id_rpcor = netcdf.defVar(ncid,'relative_power_correction','nc_float',[]);
    netcdf.putAtt(ncid,id_rpcor,'comment','flag indicating, if relative power correction (dual. pol. radars) has been used (0=yes, 1=no)');
    netcdf.putAtt(ncid,id_rpcor,'rpg_name','RelPowCorr'); % data.RelPowCorr
    
end



%%%%%%% range variables     


if isfield(data, 'Fr')
    id_Fr = netcdf.defVar(ncid,'Fr','nc_int',did_range);
    netcdf.putAtt(ncid,id_Fr,'long_name','Range factors');
    netcdf.defVarFill(ncid,id_Fr,false,int32(-999))
    netcdf.putAtt(ncid,id_Fr,'comment','See manual for details');
    netcdf.putAtt(ncid,id_Fr,'rpg_name','Fr');
end

%%%%%%%% chirp_seq_dependent variables


id_range_offsets = defh.range_offsets(ncid,did_no_seq);
netcdf.putAtt(ncid,id_range_offsets,'rpg_name','RngOffs');


id_SeqAvg = netcdf.defVar(ncid,'num_avg_chirps','nc_int',did_no_seq);
netcdf.putAtt(ncid,id_SeqAvg,'long_name','Number of averaged chirps in each chirp sequence');
netcdf.defVarFill(ncid,id_SeqAvg,false,int32(-999))
netcdf.putAtt(ncid,id_SeqAvg,'rpg_name','ChirpReps');

if isfield(data, 'ChanBW')
    id_bw = netcdf.defVar(ncid,'channel_bandwidth','nc_float',did_no_seq);
    netcdf.putAtt(ncid,id_bw,'long_name','bandwidth of individual radar channels in the sequence');
    netcdf.putAtt(ncid,id_bw,'units','Hz');
    netcdf.defVarFill(ncid,id_bw,false,NaN('single'))
    netcdf.putAtt(ncid,id_bw,'rpg_name','ChanBW');
end

if isfield(data, 'ChirpLowIF')
    id_lowIF = netcdf.defVar(ncid,'chirp_low_if','nc_float',did_no_seq);
    netcdf.putAtt(ncid,id_lowIF,'long_name','lowest IF frequency in the sequence');
    netcdf.putAtt(ncid,id_lowIF,'units','Hz');
    netcdf.defVarFill(ncid,id_lowIF,false,NaN('single'))
    netcdf.putAtt(ncid,id_lowIF,'rpg_name','ChirpLowIF');
end

if isfield(data, 'ChirpHighIF')
    id_highIF = netcdf.defVar(ncid,'chirp_high_if','nc_float',did_no_seq);
    netcdf.putAtt(ncid,id_highIF,'long_name','highest IF frequency in the sequence');
    netcdf.putAtt(ncid,id_highIF,'units','Hz');
    netcdf.defVarFill(ncid,id_highIF,false,NaN('single'))
    netcdf.putAtt(ncid,id_highIF,'rpg_name','ChirpHighIF');
end

if isfield(data, 'RangeMin')
    id_crmin = netcdf.defVar(ncid,'chirp_minrange','nc_float',did_no_seq);
    netcdf.putAtt(ncid,id_crmin,'long_name','minimum range in the sequence');
    netcdf.putAtt(ncid,id_crmin,'units','m');
    netcdf.defVarFill(ncid,id_crmin,false,NaN('single'))
    netcdf.putAtt(ncid,id_crmin,'rpg_name','RangeMin');
end

if isfield(data, 'RangeMax')
    id_crmax = netcdf.defVar(ncid,'chirp_maxrange','nc_float',did_no_seq);
    netcdf.putAtt(ncid,id_crmax,'long_name','maximum range in the sequence');
    netcdf.putAtt(ncid,id_crmax,'units','m');
    netcdf.defVarFill(ncid,id_crmax,false,NaN('single'))
    netcdf.putAtt(ncid,id_crmax,'rpg_name','RangeMax');
end

if isfield(data, 'ChirpFFTSize')
    id_rfft = netcdf.defVar(ncid,'ranging_fft','nc_float',did_no_seq);
    netcdf.putAtt(ncid,id_rfft,'long_name','Ranging FFT');
    netcdf.defVarFill(ncid,id_rfft,false,NaN('single'))
    netcdf.putAtt(ncid,id_rfft,'rpg_name','ChirpFFTSize'); % data.ChirpFFTSize
end

if isfield(data, 'ChirpInvSamples')
    id_invsmp = netcdf.defVar(ncid,'chirp_invalid_samples','nc_float',did_no_seq);
    netcdf.putAtt(ncid,id_invsmp,'long_name','number of invalid samples at beginning of chirp');
    netcdf.defVarFill(ncid,id_invsmp,false,NaN('single'))
    netcdf.putAtt(ncid,id_invsmp,'rpg_name','ChirpInvSamples'); % data.ChirpInvSamples
end

if isfield(data, 'ChirpCenterFr')
    id_ccf = netcdf.defVar(ncid,'chirp_center_freq','nc_float',did_no_seq);
    netcdf.putAtt(ncid,id_ccf,'long_name','chirp centre frequency at radar transmitter output');
    netcdf.putAtt(ncid,id_ccf,'units','MHz');
    netcdf.defVarFill(ncid,id_ccf,false,NaN('single'))
    netcdf.putAtt(ncid,id_ccf,'rpg_name','ChirpCenterFr'); % data.ChirpCenterFr
end

if isfield(data, 'ChirpBWFr')
    id_cbw = netcdf.defVar(ncid,'chirp_bandwidth','nc_float',did_no_seq);
    netcdf.putAtt(ncid,id_cbw,'long_name','chirp bandwidth at radar transmitter output');
    netcdf.putAtt(ncid,id_cbw,'units','MHz');
    netcdf.defVarFill(ncid,id_cbw,false,NaN('single'))
    netcdf.putAtt(ncid,id_cbw,'rpg_name','ChirpBWFr'); % data.ChirpBWFr
end

if isfield(data, 'FFTStartInd')
    id_ffts = netcdf.defVar(ncid,'fft_start_index','nc_float',did_no_seq);
    netcdf.putAtt(ncid,id_ffts,'long_name','start index of sequence in FFT array');
    netcdf.defVarFill(ncid,id_ffts,false,NaN('single'))
    netcdf.putAtt(ncid,id_ffts,'rpg_name','FFTStartInd'); % data.FFTStartInd
end

if isfield(data, 'FFTStopInd')
    id_ffte = netcdf.defVar(ncid,'fft_stop_index','nc_float',did_no_seq);
    netcdf.putAtt(ncid,id_ffte,'long_name','stop index of sequence in FFT array');
    netcdf.defVarFill(ncid,id_ffte,false,NaN('single'))
    netcdf.putAtt(ncid,id_ffte,'rpg_name','FFTStopInd'); % data.FFTStopInd
end

if isfield(data, 'ChirpFFTNo')
    id_fftno = netcdf.defVar(ncid,'chirp_num_fft','nc_float',did_no_seq);
    netcdf.putAtt(ncid,id_fftno,'long_name','number of FFT range layers in one chirp');
    netcdf.defVarFill(ncid,id_fftno,false,NaN('single'))
    netcdf.putAtt(ncid,id_fftno,'rpg_name','ChirpFFTNo'); % data.ChirpFFTNo
    netcdf.putAtt(ncid,id_fftno,'comment','usually 1, according to the manual');
end

id_ChrpIntTime = netcdf.defVar(ncid,'chirp_seq_int_time','nc_float',did_no_seq); % added by RG 4.12.2020
netcdf.putAtt(ncid,id_ChrpIntTime,'long_name','Integration time of each chirp sequence');
netcdf.putAtt(ncid,id_ChrpIntTime,'units','sec');
netcdf.defVarFill(ncid,id_ChrpIntTime,false,NaN('single'))
netcdf.putAtt(ncid,id_ChrpIntTime,'comment',['Calculated from chirp repetition frequency ' ...
                                             'and number of repetitions. This variable corresponds to the chirp ' ...
                                             'integration time set in the chirp table']);
                                         
id_dr = netcdf.defVar(ncid,'range_resolution','nc_float',did_no_seq);
netcdf.putAtt(ncid,id_dr,'long_name','range resolution');
netcdf.putAtt(ncid,id_dr,'units','m');
netcdf.defVarFill(ncid,id_dr,false,NaN('single'))
netcdf.putAtt(ncid,id_dr,'rpg_name','dR');


%%%%%%%% time dependend variables


id_RR = netcdf.defVar(ncid,'rain_rate','nc_float',did_time);
netcdf.putAtt(ncid,id_RR,'long_name','rain rate');
netcdf.putAtt(ncid,id_RR,'units','mm h-1');
netcdf.defVarFill(ncid,id_RR,false,NaN('single'))
netcdf.putAtt(ncid,id_RR,'comment','measured by the meteo-station attached to the radar');
netcdf.putAtt(ncid,id_RR,'rpg_name','RR');

id_rh = netcdf.defVar(ncid,'rh','nc_float',did_time);
netcdf.putAtt(ncid,id_rh,'long_name','relative humidity');
netcdf.putAtt(ncid,id_rh,'units','%');
netcdf.defVarFill(ncid,id_rh,false,NaN('single'))
netcdf.putAtt(ncid,id_rh,'comment','measured by the meteo-station attached to the radar');
netcdf.putAtt(ncid,id_rh,'rpg_name','rh');   

id_T_env = netcdf.defVar(ncid,'surface_temperature','nc_float',did_time);
netcdf.putAtt(ncid,id_T_env,'long_name','ambient air temperature');
netcdf.putAtt(ncid,id_T_env,'units','K');
netcdf.defVarFill(ncid,id_T_env,false,NaN('single'))
netcdf.putAtt(ncid,id_T_env,'comment','measured by the meteo-station attached to the radar');
netcdf.putAtt(ncid,id_T_env,'rpg_name','T_env');

id_pres = netcdf.defVar(ncid,'pressure','nc_float',did_time);
netcdf.putAtt(ncid,id_pres,'long_name','surface pressure');
netcdf.putAtt(ncid,id_pres,'units','hPa');
netcdf.defVarFill(ncid,id_pres,false,NaN('single'))
netcdf.putAtt(ncid,id_pres,'comment','measured by the meteo-station attached to the radar');
netcdf.putAtt(ncid,id_pres,'rpg_name','press');

id_ff = netcdf.defVar(ncid,'wspeed','nc_float',did_time);
netcdf.putAtt(ncid,id_ff,'long_name','wind speed');
netcdf.putAtt(ncid,id_ff,'units','m s-1');
netcdf.defVarFill(ncid,id_ff,false,NaN('single'))
netcdf.putAtt(ncid,id_ff,'comment','measured by the meteo-station attached to the radar');
netcdf.putAtt(ncid,id_ff,'rpg_name','ff');

id_fff = netcdf.defVar(ncid,'wdir','nc_float',did_time);
netcdf.putAtt(ncid,id_fff,'long_name','wind direction');
netcdf.putAtt(ncid,id_fff,'units','degrees');
netcdf.defVarFill(ncid,id_fff,false,NaN('single'))
netcdf.putAtt(ncid,id_fff,'comment','wind direction follows standard meteorological convention (0 degrees for wind from North); measured by the meteo-station attached to the radar');
netcdf.putAtt(ncid,id_fff,'rpg_name','fff');

id_vol = netcdf.defVar(ncid,'volt','nc_float',did_time);
netcdf.putAtt(ncid,id_vol,'long_name','Voltage direct detection channel');
netcdf.putAtt(ncid,id_vol,'units','V');
netcdf.defVarFill(ncid,id_vol,false,NaN('single'))
netcdf.putAtt(ncid,id_vol,'rpg_name','DDVolt');

id_Tb = netcdf.defVar(ncid,'tb','nc_float',did_time);
netcdf.putAtt(ncid,id_Tb,'long_name','brightness temperature 89 GHz');
netcdf.putAtt(ncid,id_Tb,'units','K');
netcdf.defVarFill(ncid,id_Tb,false,NaN('single'))
netcdf.putAtt(ncid,id_Tb,'comment','measurement of the 89-GHz channel');
netcdf.putAtt(ncid,id_Tb,'rpg_name','Tb');

id_lwp = netcdf.defVar(ncid,'lwp','nc_float',did_time);
netcdf.putAtt(ncid,id_lwp,'long_name','liquid water path');
netcdf.putAtt(ncid,id_lwp,'units','g m-2');
netcdf.defVarFill(ncid,id_lwp,false,NaN('single'))
netcdf.putAtt(ncid,id_lwp,'comment',['The liquid water path is calculated from '...
                                     'the tb measurement of the 89-GHz channel and '...
                                     'provided by the instrument software using a retrieval '...                                     
                                     'that was developed by RPG and '...
                                     'is based on a neural network approach']);
netcdf.putAtt(ncid,id_lwp,'rpg_name','lwp');


id_powIF = netcdf.defVar(ncid,'pow_if','nc_float',did_time);
netcdf.putAtt(ncid,id_powIF,'long_name','IF power at ADC');
netcdf.putAtt(ncid,id_powIF,'units','microWatt');
netcdf.defVarFill(ncid,id_powIF,false,NaN('single'))
netcdf.putAtt(ncid,id_powIF,'rpg_name','PowIF');

id_ele = netcdf.defVar(ncid,'scanner_elevation','nc_float',did_time);
netcdf.putAtt(ncid,id_ele,'long_name','Scanner elevation angle');
netcdf.putAtt(ncid,id_ele,'units','degrees');
netcdf.defVarFill(ncid,id_ele,false,NaN('single'))
netcdf.putAtt(ncid,id_ele,'rpg_name','Elev'); % "elevation angle of sample"

id_az = netcdf.defVar(ncid,'scanner_azimuth','nc_float',did_time);
netcdf.putAtt(ncid,id_az,'long_name','Scanner azimuth angle');
netcdf.putAtt(ncid,id_az,'units','degrees');
netcdf.defVarFill(ncid,id_az,false,NaN('single'))
netcdf.putAtt(ncid,id_az,'rpg_name','Azi');  % "aszimuth angle of sample"

                                                                 
id_T_trans = netcdf.defVar(ncid,'temp_trans','nc_float',did_time);
netcdf.putAtt(ncid,id_T_trans,'long_name','transmitter temperature');
netcdf.putAtt(ncid,id_T_trans,'units','K');
netcdf.defVarFill(ncid,id_T_trans,false,NaN('single'))
netcdf.putAtt(ncid,id_T_trans,'comment',['For a stable performance 309K < t_trans < 313K.']);
netcdf.putAtt(ncid,id_T_trans,'rpg_name','T_trans');

                                 
id_T_rec = netcdf.defVar(ncid,'temp_rec','nc_float',did_time);
netcdf.putAtt(ncid,id_T_rec,'long_name','receiver temperature');
netcdf.putAtt(ncid,id_T_rec,'units','K');
netcdf.defVarFill(ncid,id_T_rec,false,NaN('single'))
netcdf.putAtt(ncid,id_T_rec,'comment',['For a stable performance 306K < t_rec < 312K.']);
netcdf.putAtt(ncid,id_T_rec,'rpg_name','T_rec');


id_T_pc = netcdf.defVar(ncid,'temp_pc','nc_float',did_time);
netcdf.putAtt(ncid,id_T_pc,'long_name','pc temperature');
netcdf.putAtt(ncid,id_T_pc,'units','K');
netcdf.defVarFill(ncid,id_T_pc,false,NaN('single'))
netcdf.putAtt(ncid,id_T_pc,'comment',['For a stable performance t_pc < 323K.']);
netcdf.putAtt(ncid,id_T_pc,'rpg_name','T_pc');

if isfield(data, 'reserved')
    id_incel = netcdf.defVar(ncid,'instrument_elevation','nc_float',did_time); 
    netcdf.putAtt(ncid,id_incel,'long_name','instrument elevation angle');
    netcdf.putAtt(ncid,id_incel,'units','degrees');
    netcdf.defVarFill(ncid,id_incel,false,NaN('single'))
    netcdf.putAtt(ncid,id_incel,'comment','recorded by an internal inclination sensor, if included in the instrument');

    id_incea = netcdf.defVar(ncid,'instrument_azimuth','nc_float',did_time);
    netcdf.putAtt(ncid,id_incea,'long_name','instrument azimuth angle');
    netcdf.putAtt(ncid,id_incea,'units','degrees');
    netcdf.defVarFill(ncid,id_incea,false,NaN('single'))
    netcdf.putAtt(ncid,id_incea,'comment','recorded by an internal inclination sensor, if included in the instrument');
end


if isfield(data, 'QF')
    id_QF = netcdf.defVar(ncid,'rpg_quality_flag','nc_byte',did_time); 
    netcdf.putAtt(ncid,id_QF,'long_name','quality flag from RPG software');
    netcdf.putAtt(ncid,id_QF,'comment','Bit 1: ADC saturation, Bit 2: spectral width too high, Bit 3: no transm. power leveling');
    netcdf.putAtt(ncid,id_QF,'rpg_name','QF');
end


id_timeshift = netcdf.defVar(ncid,'timestamp_correction','nc_float',did_time); 
netcdf.putAtt(ncid,id_timeshift,'long_name','time stamp correction');
netcdf.putAtt(ncid,id_timeshift,'units','sec');
netcdf.defVarFill(ncid,id_timeshift,false,NaN('single'))
netcdf.putAtt(ncid,id_timeshift,'comment','Sanity checks correct the time stamp if duplicate time stamps or backward moving time detected. This variable contains the magnitude of the time stamp correction. Get the original time stamp by calculating (time+sampleTms*1e-3) - timestamp_correction');

                                         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% multi-D variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


id_Aliasmask = netcdf.defVar(ncid,'alias_mask','nc_byte',[did_range,did_time]);
netcdf.putAtt(ncid,id_Aliasmask,'long_name','aliasing mask');
netcdf.putAtt(ncid,id_Aliasmask,'comment',['AliasMask indicates the bins where aliasing is detected. '...
                                           '0 = no aliasing; 1 = aliasing occurs. Only applicable if '...
                                           'variable alias_method_flag is 2.']);

id_AliasStatus = netcdf.defVar(ncid,'dealias_status','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_AliasStatus,'long_name','Flag indicating quality of dealiasing.');
netcdf.putAtt(ncid,id_AliasStatus,'comment',['Must be converted into four bit binary string. If 0, i.e. dec2bin(AliasStatus) = 0000 then dealiasing was successful. ',...
    'If 2^0 bit is 1, there was no initial guess velocity for this bin. ',...
    'If 2^1 bit is 1, the spectrum where the main peak was assumed to occur could not be concatenated properly since chirp sequence boundaries were reached. ',...
    'If 2^2 bit is 1, there is still significant signal close to the nyquist limit. dealiasing might not have been performed properly. ',...
    'If 2^3 bit is 1, the mean of the mean Doppler velocity in differs by more than 5 m/s from the neighbouring bins mean value. ',...
    'The larger the decimal number, the less reliable the data. ',...
    'This flag is only available if dealiasing was applied in process_joyrad94_data.m. For RPG dealising there is no quality indicator.']);


id_MinVel_Correction = netcdf.defVar(ncid,'MinVel_Correction','nc_float',[did_range,did_time]);
netcdf.putAtt(ncid,id_MinVel_Correction,'long_name','Minimum velocity correction');
netcdf.putAtt(ncid,id_MinVel_Correction,'units','m s-1');
netcdf.defVarFill(ncid,id_MinVel_Correction,false,NaN('single'))
netcdf.putAtt(ncid,id_MinVel_Correction,'comment',['If dealising was applied the array indicates by how much MinVel '...
    'was corrected by the final quality check (output from dealias_spectra_vm_column_quality_check.m). '...
    'Subtracting this value from MinVel provides the offset before the final quality check (output from dealias_spectra.m).']);

if isfield(data, 'PNv')
    id_PNv = netcdf.defVar(ncid,'PNv','nc_float',[did_range,did_time]);
    netcdf.putAtt(ncid,id_PNv,'long_name','Total IF power in vertical polarization measured at ADC input');
    netcdf.defVarFill(ncid,id_PNv,false,NaN('single'))
    netcdf.putAtt(ncid,id_PNv,'rpg_name','PNv');
end

if isfield(data, 'lineart_flag')
    id_lineartf = netcdf.defVar(ncid,'flag_artifact','nc_byte',did_range);
    netcdf.putAtt(ncid,id_lineartf,'long_name','flag for so-called line artifact');
    netcdf.putAtt(ncid,id_lineartf,'flag_values', [0b0, 0b1]);
    netcdf.putAtt(ncid,id_lineartf,'flag_meanings','data_intact data_affected');
    netcdf.putAtt(ncid,id_lineartf,'comment','variable indicates the range gates where data is affected, and partially corrected');
    netcdf.defVarFill(ncid,id_lineartf,false,0b10)
end

if Temp_ret_flag
    id_tprof = netcdf.defVar(ncid,'temperature_profile','nc_float',[did_T_range,did_time]);
    netcdf.putAtt(ncid,id_tprof,'long_name','temperature profile');
    netcdf.defVarFill(ncid,id_tprof,false,NaN('single'))
    netcdf.putAtt(ncid,id_tprof,'comment','retrieved by RPG software');
    netcdf.putAtt(ncid,id_tprof,'rpg_name','TPr');
end

if Hum_ret_flag
    id_hprof = netcdf.defVar(ncid,'abshum_profile','nc_float',[did_H_range,did_time]);
    netcdf.putAtt(ncid,id_hprof,'long_name','absolute humidity profile');
    netcdf.defVarFill(ncid,id_hprof,false,NaN('single'))
    netcdf.putAtt(ncid,id_hprof,'comment','retrieved by RPG software');
    netcdf.putAtt(ncid,id_hprof,'rpg_name','AHPr');

    id_rhprof = netcdf.defVar(ncid,'RH_profile','nc_float',[did_H_range,did_time]);
    netcdf.putAtt(ncid,id_rhprof,'long_name','relative humidity profile');
    netcdf.defVarFill(ncid,id_rhprof,false,NaN('single'))
    netcdf.putAtt(ncid,id_rhprof,'comment','retrieved by RPG software');
    netcdf.putAtt(ncid,id_rhprof,'rpg_name','RHPr');
end

if data.DualPol > 0
    disp('WARMING!!! No polarimetric variables included in the output files')
end


%% ###################### compression

% range dependent variables
netcdf.defVarDeflate(ncid,id_height,true,true,9);

if isfield(data, 'Fr')
    netcdf.defVarDeflate(ncid,id_Fr,true,true,9);
end

% chirp_seq dependent variables
netcdf.defVarDeflate(ncid,id_range_offsets,true,true,9);
netcdf.defVarDeflate(ncid,id_SeqAvg,true,true,9);
if isfield(data, 'ChanBW')
    netcdf.defVarDeflate(ncid,id_bw,true,true,9);
end
if isfield(data, 'ChirpLowIF')
    netcdf.defVarDeflate(ncid,id_lowIF,true,true,9);
end
if isfield(data, 'ChirpHighIF')
    netcdf.defVarDeflate(ncid,id_highIF,true,true,9);
end
if isfield(data, 'RangeMin')
    netcdf.defVarDeflate(ncid,id_crmin,true,true,9);
end
if isfield(data, 'RangeMax')
    netcdf.defVarDeflate(ncid,id_crmax,true,true,9);
end
if isfield(data, 'ChirpFFTSize')
    netcdf.defVarDeflate(ncid,id_rfft,true,true,9);
end
if isfield(data, 'ChirpInvSamples')
    netcdf.defVarDeflate(ncid,id_invsmp,true,true,9);
end
if isfield(data, 'ChirpCenterFr')
    netcdf.defVarDeflate(ncid,id_ccf,true,true,9);
end
if isfield(data, 'ChirpBWFr')
    netcdf.defVarDeflate(ncid,id_cbw,true,true,9);
end
if isfield(data, 'FFTStartInd')
    netcdf.defVarDeflate(ncid,id_ffts,true,true,9);
end
if isfield(data, 'FFTStopInd')
    netcdf.defVarDeflate(ncid,id_ffte,true,true,9);
end
if isfield(data, 'ChirpFFTNo')
    netcdf.defVarDeflate(ncid,id_fftno,true,true,9);
end
netcdf.defVarDeflate(ncid,id_ChrpIntTime,true,true,9);
netcdf.defVarDeflate(ncid,id_dr,true,true,9);

% time dependend variables
netcdf.defVarDeflate(ncid,id_sampleTms,true,true,9);
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
netcdf.defVarDeflate(ncid,id_T_trans,true,true,9);
netcdf.defVarDeflate(ncid,id_T_rec,true,true,9);
netcdf.defVarDeflate(ncid,id_T_pc,true,true,9);
if isfield(data, 'reserved')
    netcdf.defVarDeflate(ncid,id_incel,true,true,9);
    netcdf.defVarDeflate(ncid,id_incea,true,true,9);
end
if isfield(data, 'QF')
    netcdf.defVarDeflate(ncid,id_QF,true,true,9);
end

netcdf.defVarDeflate(ncid,id_timeshift,true,true,9);


% multidimensional variables
netcdf.defVarDeflate(ncid,id_Aliasmask,true,true,9);
netcdf.defVarDeflate(ncid,id_AliasStatus,true,true,9);
netcdf.defVarDeflate(ncid,id_MinVel_Correction,true,true,9);
if isfield(data, 'PNv')
    netcdf.defVarDeflate(ncid,id_PNv,true,true,9);
end
if isfield(data, 'lineart_flag')
    netcdf.defVarDeflate(ncid,id_lineartf,true,true,9);
end

if Temp_ret_flag
    netcdf.defVarDeflate(ncid,id_tprof,true,true,9);
end
if Hum_ret_flag
    netcdf.defVarDeflate(ncid,id_hprof,true,true,9);
    netcdf.defVarDeflate(ncid,id_rhprof,true,true,9);
end

netcdf.endDef(ncid);

%% ####################### put variables into file

% variables for dimensions
netcdf.putVar(ncid,id_time,0,data.totsamp,data.time);
netcdf.putVar(ncid,id_range,0,data.n_levels,data.range);
netcdf.putVar(ncid,id_no_seq,0,data.no_chirp_seq,1:data.no_chirp_seq);

if Temp_ret_flag
    netcdf.putVar(ncid,id_coord_Tlevels,0,data.T_altcount,data.T_alt);
end
if Hum_ret_flag
    netcdf.putVar(ncid,id_coord_Hlevels,0,data.H_altcount,data.H_alt);
end

% scalars
netcdf.putVar(ncid,id_filecode,data.filecode);
netcdf.putVar(ncid,id_progno,data.progno);
if isfield(data, 'custname')
    netcdf.putVar(ncid,id_custname,0,length(data.custname), data.custname);
end
% netcdf.putVar(ncid,id_samples,data.totsamp);
% netcdf.putVar(ncid,id_levels,data.n_levels);
% netcdf.putVar(ncid,id_no_chirp_seq,data.no_chirp_seq);
netcdf.putVar(ncid,id_AntiAlias,data.AntiAlias);
netcdf.putVar(ncid,id_pol,data.DualPol);

if Temp_ret_flag
    netcdf.putVar(ncid,id_T_levels,data.T_altcount);
end
if Hum_ret_flag
    netcdf.putVar(ncid,id_H_levels,data.H_altcount);
end
netcdf.putVar(ncid,id_compress,data.CompEna);
% netcdf.putVar(ncid,id_cal_mom,data.cal_mom);
netcdf.putVar(ncid,id_swv,0,length(data.radarsw), data.radarsw);
netcdf.putVar(ncid,id_AntSep,data.AntSep);
netcdf.putVar(ncid,id_AntDia,data.AntDia);

netcdf.putVar(ncid,id_C,data.C);
if isfield(data, 'SampRate')
    netcdf.putVar(ncid,id_smprt,data.SampRate);
end
if isfield(data, 'MaxRange')
    netcdf.putVar(ncid,id_maxr,data.MaxRange);
end
if isfield(data, 'SupPowLev')
    netcdf.putVar(ncid,id_powesup,data.SupPowLev);
end
if isfield(data, 'SpkFilEna')
    netcdf.putVar(ncid,id_sfil,data.SpkFilEna);
end
if isfield(data, 'FFTWindow')
    FFTWindows = {'SQUARE', 'PARZEN', 'BLACKMAN', 'WELCH', 'SLEPIAN2', 'SLEPIAN3'};
    netcdf.putVar(ncid,id_fftwin,0,length(FFTWindows{data.FFTWindow+1}),FFTWindows{data.FFTWindow+1});
    % 0 = SQUARE
    % 1 = PARZEN
    % 2 = BLACKMAN
    % 3 = WELCH
    % 4 = SLEPIAN2
    % 5 = SLEPIAN3
end
if isfield(data, 'FFTInputRng')
    netcdf.putVar(ncid,id_finrn,data.FFTInputRng);
end
if isfield(data, 'NoiseFilt')
    netcdf.putVar(ncid,id_noisfil,data.NoiseFilt);
end
netcdf.putVar(ncid,id_SampDur,data.SampDur);

% range dependent
netcdf.putVar(ncid,id_height,0,data.n_levels,data.range+config.MSL);
if isfield(data, 'Fr')
    netcdf.putVar(ncid,id_Fr,0,data.n_levels,data.Fr);
end

% chrip seq dependent variables
netcdf.putVar(ncid,id_range_offsets,0,data.no_chirp_seq,data.range_offsets-1);
netcdf.putVar(ncid,id_SeqAvg,0,data.no_chirp_seq,data.SeqAvg);
if isfield(data, 'ChanBW')
    netcdf.putVar(ncid,id_bw,0,data.no_chirp_seq,data.ChanBW);
end
if isfield(data, 'ChirpLowIF')
    netcdf.putVar(ncid,id_lowIF,0,data.no_chirp_seq,data.ChirpLowIF);
end
if isfield(data, 'ChirpHighIF')
    netcdf.putVar(ncid,id_highIF,0,data.no_chirp_seq,data.ChirpHighIF);
end
if isfield(data, 'RangeMin')
    netcdf.putVar(ncid,id_crmin,0,data.no_chirp_seq,data.RangeMin);
end
if isfield(data, 'RangeMax')
    netcdf.putVar(ncid,id_crmax,0,data.no_chirp_seq,data.RangeMax);
end
if isfield(data, 'ChirpFFTSize')
    netcdf.putVar(ncid,id_rfft,0,data.no_chirp_seq,data.ChirpFFTSize);
end
if isfield(data, 'ChirpInvSamples')
    netcdf.putVar(ncid,id_invsmp,0,data.no_chirp_seq,data.ChirpInvSamples);
end
if isfield(data, 'ChirpCenterFr')
    netcdf.putVar(ncid,id_ccf,0,data.no_chirp_seq,data.ChirpCenterFr);
end
if isfield(data, 'ChirpBWFr')
    netcdf.putVar(ncid,id_cbw,0,data.no_chirp_seq,data.ChirpBWFr);
end
if isfield(data, 'FFTStartInd')
    netcdf.putVar(ncid,id_ffts,0,data.no_chirp_seq,data.FFTStartInd);
end
if isfield(data, 'FFTStopInd')
    netcdf.putVar(ncid,id_ffte,0,data.no_chirp_seq,data.FFTStopInd)
end
if isfield(data, 'ChirpFFTNo')
    netcdf.putVar(ncid,id_fftno,0,data.no_chirp_seq,data.ChirpFFTNo)
end

netcdf.putVar(ncid,id_ChrpIntTime,0,data.no_chirp_seq,data.ChirpIntTime); % added by RG 4.12.2020
netcdf.putVar(ncid,id_dr,0,data.no_chirp_seq,data.dr);

% time dependent variables

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
netcdf.putVar(ncid,id_T_trans,0,data.totsamp,data.T_trans);
netcdf.putVar(ncid,id_T_rec,0,data.totsamp,data.T_rec);
netcdf.putVar(ncid,id_T_pc,0,data.totsamp,data.T_pc);
if isfield(data, 'reserved')
    netcdf.putVar(ncid,id_incel,0,data.totsamp,data.reserved(:,2));
    netcdf.putVar(ncid,id_incea,0,data.totsamp,data.reserved(:,3));
end
if isfield(data, 'QF')
    netcdf.putVar(ncid,id_QF,0,data.totsamp,data.QF);
end

netcdf.putVar(ncid,id_timeshift,0,data.totsamp,data.timeshift);

% multidimensional variables
netcdf.putVar(ncid,id_Aliasmask,[0,0],[data.n_levels,data.totsamp],data.Aliasmask');
netcdf.putVar(ncid,id_AliasStatus,[0,0],[data.n_levels,data.totsamp],data.AliasStatus');
netcdf.putVar(ncid,id_MinVel_Correction,[0,0],[data.n_levels,data.totsamp],data.MinVel_Correction');
if isfield(data, 'PNv')
    netcdf.putVar(ncid,id_PNv,[0,0],[data.n_levels,data.totsamp],data.PNv');
end
if isfield(data, 'lineart_flag')
    netcdf.putVar(ncid,id_lineartf,0,data.n_levels, double(data.lineart_flag) );
end
if Temp_ret_flag
    netcdf.putVar(ncid,id_tprof,[0,0],[data.T_altcount, data.totsamp],data.Tprof');
end
if Hum_ret_flag
    netcdf.putVar(ncid,id_hprof,[0,0],[data.H_altcount,data.totsamp],data.Qprof');
    netcdf.putVar(ncid,id_rhprof,[0,0],[data.H_altcount,data.totsamp],data.RHprof');
end

netcdf.close(ncid);


end % function
