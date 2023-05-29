% collection of functions to make it easy to use same attributes for output
% variables in different netcdf files
% RG 22.9.2022

function fh = outvarmeta

    fh.height = @height;
    fh.msl = @msl;
    fh.range = @range;
    fh.time = @time;
    fh.lat = @lat;
    fh.lon = @lon;
    fh.sampleTms = @sampleTms;
    fh.chirp_sequence = @chirp_sequence;
    fh.radar_software = @radar_software;
    fh.radar_software_as_str = @radar_software_as_str;
    fh.DoppMax = @DoppMax;
    fh.SeqIntTime = @SeqIntTime;
    fh.range_offsets = @range_offsets;
    fh.SLv = @SLv;
    fh.noisestd = @noisestd; 
    fh.zecalib = @zecalib; 
    fh.aggregFlag = @aggregFlag; 
    fh.ze_comment = @ze_comment;
    
end % fh



function id_height = height(ncid, did_height)

    id_height = netcdf.defVar(ncid,'height','nc_float',did_height);
    netcdf.putAtt(ncid,id_height,'long_name',['height above mean sea level, at the '...
                                            'center of the radar range gate']);
    netcdf.putAtt(ncid,id_height,'units','m');
    netcdf.putAtt(ncid,id_height,'positive','up'); % recommended by CF convention
    netcdf.putAtt(ncid,id_height,'comment','calculated as range + instrument_altitude');
    
end


function id_MSL = msl(ncid)

    id_MSL = netcdf.defVar(ncid,'instrument_altitude','nc_float',[]);
    netcdf.putAtt(ncid,id_MSL,'long_name','instrument altitude above mean sea level');
    netcdf.putAtt(ncid,id_MSL,'units','m');
    netcdf.defVarFill(ncid,id_MSL,false,NaN('single'))
    
end



function id_range = range(ncid, did_range)
    id_range = netcdf.defVar(ncid,'range','nc_float',did_range);
    netcdf.putAtt(ncid,id_range,'long_name','range');
    netcdf.putAtt(ncid,id_range,'units','m');
    netcdf.putAtt(ncid,id_range,'comment',['Range from the radar antenna to the '...
                                           'center of the radar range gate']);
end
                                   

function id_time = time(ncid, did_time, timetype, totsamp_flag)

    id_time = netcdf.defVar(ncid,'time',timetype,did_time);
    netcdf.putAtt(ncid,id_time,'standard_name','time');
    netcdf.putAtt(ncid,id_time,'units','seconds since 2001-01-01');
    netcdf.putAtt(ncid,id_time,'calendar','standard');
    netcdf.putAtt(ncid,id_time,'comment','timestamp given at the end of the sample');
    
    if strcmp(timetype, 'nc_int') || strcmp(timetype, 'NC_INT')
        netcdf.putAtt(ncid,id_time,'ancillary_variables','sample_tms')
    end
    if totsamp_flag
        netcdf.putAtt(ncid,id_time, 'note', ['Issues with timestamps occurred (value outside the expected range, '...
                                             'time moving backwards, and/or dublicate timestamps found in the lv0-file). ' ...
                                             'Depending on the issue, profiles removed or timestamp edited.'])
    end
    
end


function id_lat = lat(ncid)
    id_lat = netcdf.defVar(ncid,'lat','nc_float',[]);
    netcdf.putAtt(ncid,id_lat,'long_name','latitude');
    netcdf.putAtt(ncid,id_lat,'standard_name','latitude');
    netcdf.putAtt(ncid,id_lat,'units','degrees_north');
    netcdf.defVarFill(ncid,id_lat,false,NaN('single'))
end


function id_lon = lon(ncid)
    id_lon = netcdf.defVar(ncid,'lon','nc_float',[]);
    netcdf.putAtt(ncid,id_lon,'long_name','longitude');
    netcdf.putAtt(ncid,id_lon,'standard_name','longitude');
    netcdf.putAtt(ncid,id_lon,'units','degrees_east');
    netcdf.defVarFill(ncid,id_lon,false,NaN('single'))
end

function id_sampleTms = sampleTms(ncid, did_time)
    id_sampleTms = netcdf.defVar(ncid,'sample_tms','nc_int',did_time);
    netcdf.putAtt(ncid,id_sampleTms,'long_name','milliseconds of sample');
    netcdf.putAtt(ncid,id_sampleTms,'units','millisec');
    netcdf.defVarFill(ncid,id_sampleTms,false,int32(-999))
    netcdf.putAtt(ncid,id_sampleTms,'comment','add to time to get the precise sample time');  

end

function id_no_seq = chirp_sequence(fileid, did_no_seq)
    id_no_seq = netcdf.defVar(fileid,'chirp_sequence','nc_int',did_no_seq);
    netcdf.putAtt(fileid,id_no_seq,'long_name','chirp sequences');
    netcdf.putAtt(fileid,id_no_seq,'units','count');
%     netcdf.defVarFill(fileid,id_no_seq,false,int32(-999)) % coordinate
%     variables are not supposed to have missing values
end


function id_swv = radar_software(fileid)
    id_swv = netcdf.defVar(fileid,'radar_software','nc_float',[]);
    netcdf.putAtt(fileid,id_swv,'long_name','radar software version');
    netcdf.defVarFill(fileid,id_swv,false,NaN('single'))
end


function id_swv = radar_software_as_str(fileid, did_str)
    id_swv = netcdf.defVar(fileid,'radar_software','nc_char', did_str);
    netcdf.putAtt(fileid,id_swv,'long_name','radar software version');
end


function id_DoppMax = DoppMax(fileid, did_no_seq)
    id_DoppMax = netcdf.defVar(fileid,'nyquist_velocity','nc_float',did_no_seq);
    netcdf.putAtt(fileid,id_DoppMax,'long_name',['Nyquist velocity']);
    netcdf.putAtt(fileid,id_DoppMax,'units','m s-1');
    netcdf.defVarFill(fileid,id_DoppMax,false,NaN('single'))
    netcdf.putAtt(fileid,id_DoppMax,'comment','max. unambiguous Doppler velocity; due to dealiasing values outside of the +-Nyquist velocity range are possible');
end

function id_SeqIntTime = SeqIntTime(fileid, did_no_seq)
    id_SeqIntTime = netcdf.defVar(fileid,'chirpseq_effinttime','nc_float',did_no_seq);
    netcdf.putAtt(fileid,id_SeqIntTime,'long_name','Effective integration time of each chirp sequence');
    netcdf.putAtt(fileid,id_SeqIntTime,'units','sec');
    netcdf.defVarFill(fileid,id_SeqIntTime,false,NaN('single'))
    netcdf.putAtt(fileid,id_SeqIntTime,'comment','Provided by radar software.');
    
end

function id_range_offsets = range_offsets(fileid,did_no_seq)
    id_range_offsets = netcdf.defVar(fileid,'chirpseq_startix','nc_int',did_no_seq); % used to be range_offsets
    netcdf.putAtt(fileid,id_range_offsets,'long_name','chirp sequence start index');
    netcdf.defVarFill(fileid,id_range_offsets,false,int32(-999))
    netcdf.putAtt(fileid,id_range_offsets,'comment','index for the height layer where each chirp sequence starts; zero-based indexing');
end


function id_SLv = SLv(fileid, did_height, did_time)
    id_SLv = netcdf.defVar(fileid,'noise_threshold','nc_float',[did_height,did_time]);
    netcdf.putAtt(fileid,id_SLv,'long_name','signal strength threshold used for data logging');
    netcdf.putAtt(fileid,id_SLv,'units','mm6/m3');
    netcdf.defVarFill(fileid,id_SLv,false,NaN('single'))
    netcdf.putAtt(fileid,id_SLv,'comment',['In the radar software data is only logged '...
                                           'for a range bin where signal exceeds a '...
                                           'threshold in at least one Doppler bin at the '...
                                           'given range. The threshold is defined as a ' ...
                                           'multiple of the standard deviation of the ' ...
                                           'spectrally-resolved noise power. noise_threshold ' ...
                                           'gives the threshold calibrated in mm6/m3 as '...
                                           'provided by the radar software, and it '...
                                           'determines the sensitivity limit of the '...
                                           'measurement. The calibration correction, '...
                                           'given by variable ze_calibration, has not '...
                                           'been added to noise_threshold.']);
end

function id_NStd = noisestd(fileid, did_no_seq, did_time)
    id_NStd = netcdf.defVar(fileid,'noise_threshold','nc_float',[did_no_seq,did_time]);
    netcdf.putAtt(fileid,id_NStd,'long_name','signal strength threshold used for data logging');
    netcdf.putAtt(fileid,id_NStd,'units','mm6/m3');
    netcdf.defVarFill(fileid,id_NStd,false,NaN('single'))
    netcdf.putAtt(fileid,id_NStd,'comment','In the radar software data is only logged for a range bin where signal exceeds a threshold in at least one Doppler bin at the given range. The threshold is defined as the standard deviation of noise power in each chirp sequence. Note, that noise_threshold is calculated differently in later radar software versions. noise_threshold gives the threshold as provided by the radar software, and it determines the sensitivity limit of the measurement. The calibration correction, given by variable ze_calibration, has not been added to noise_threshold.');
end

function id_ZeCalib = zecalib(fileid)
    id_ZeCalib = netcdf.defVar(fileid,'ze_calibration','nc_float',[]);
    netcdf.putAtt(fileid,id_ZeCalib,'long_name','Reflectivity calibration correction');
    netcdf.putAtt(fileid,id_ZeCalib,'units','dB');
    netcdf.defVarFill(fileid,id_ZeCalib,false,NaN('single'))
    netcdf.putAtt(fileid,id_ZeCalib,'comment','Calibration correction added to radar reflectivity. To get uncorrected reflectivity: ze_unccorected = ze - ze_calibration. Missing value indicates no correction applied. Radar calibration is corrected using a reference radar, for further details see documentation.');
end

function id_QF = aggregFlag(fileid, did_time, RPGflag)

    id_QF = netcdf.defVar(fileid,'quality_flag','nc_byte',did_time);
    netcdf.putAtt(fileid,id_QF,'long_name','quality flag');
    netcdf.putAtt(fileid,id_QF,'standard_name','aggregate_quality_flag');
    
    if RPGflag
        netcdf.putAtt(fileid,id_QF,'flag_masks',[0b0001, 0b0010, 0b0100, 0b1000]);
        netcdf.putAtt(fileid,id_QF,'flag_meanings','timestamp_flag dealiasing_flag dealias_problem_flag radar_operation_flag');
        netcdf.putAtt(fileid,id_QF,'comment',['Aggregate flag to collect the most important information on data quality. ' ...
                                            'Bit0 (2^0): timestamp has been edited in data processing. '...
                                            'Bit1 (2^1): column has been dealiased (see documentation for details). '...
                                            'Bit2 (2^2): issues in dealiasing, larger uncertainty in vm. '...
                                            'Bit3 (2^3): any of the flags related to radar operation set in radar software (ADC saturation, spectral width too high, no transmitter power leveling), measurement quality possibly impacted.'...
                                            ]);
    
    else  % in software version 1 no RPG flag available

        netcdf.putAtt(fileid,id_QF,'flag_masks',[0b0001, 0b0010, 0b0100]);
        netcdf.putAtt(fileid,id_QF,'flag_meanings','timestamp_flag dealiasing_flag dealias_problem_flag');
        netcdf.putAtt(fileid,id_QF,'comment',['Aggregate flag to collect the most important information on data quality. ' ...
                                            'Bit0 (2^0): timestamp has been edited in data processing. '...
                                            'Bit1 (2^1): column has been dealiased (see documentation for details). '...
                                            'Bit2 (2^2): issues in dealiasing, larger uncertainty in vm. '...
                                            ]);
    end
end


function ze_comment(data, ncid, id_Ze)

    Ze_calibcorr_label = 'If a calibration correction has been applied, it is already included in this variable, and the value is given by variable ze_calibration.';

    if isfield(data, 'Ze_label') % Ze corrected, adding note
        netcdf.putAtt(ncid,id_Ze,'comment', [data.Ze_label ' ' Ze_calibcorr_label ]);
        netcdf.putAtt(ncid,id_Ze,'corretion_dB',data.Ze_corr);
    else
        netcdf.putAtt(ncid,id_Ze,'comment', Ze_calibcorr_label);
    end

end 
