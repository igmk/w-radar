% collection of functions to make it easy to use same attributes for output
% variables in different netcdf files
% RG 22.9.2022

function fh = outvarmeta

    fh.height = @height;
    fh.range = @range;
    fh.time = @time;
    fh.lat = @lat;
    fh.lon = @lon;
    fh.sampleTms = @sampleTms;
    fh.chirp_sequence = @chirp_sequence;
    fh.radar_software = @radar_software;
    fh.DoppMax = @DoppMax;
    fh.SeqIntTime = @SeqIntTime;
    fh.range_offsets = @range_offsets;
    fh.SLv = @SLv;
    
end % fh



function id_height = height(ncid, did_height)

    id_height = netcdf.defVar(ncid,'height','nc_float',did_height);
    netcdf.putAtt(ncid,id_height,'long_name',['height above mean sea level, at the '...
                                            'center of the radar range gate']);
    netcdf.putAtt(ncid,id_height,'units','m');
    netcdf.putAtt(ncid,id_height,'positive','up'); % recommended by CF convention
    netcdf.putAtt(ncid,id_height,'comment','calculated as range + instrument_altitude');
    
end


function id_range = range(ncid, did_range)
    id_range = netcdf.defVar(ncid,'range','nc_float',did_range);
    netcdf.putAtt(ncid,id_range,'long_name','range');
    netcdf.putAtt(ncid,id_range,'units','m');
    netcdf.putAtt(ncid,id_range,'comment',['Range from the radar antenna to the '...
                                           'center of the radar range gate']);
end
                                   

function id_time = time(ncid, did_time, totsamp_flag)

    id_time = netcdf.defVar(ncid,'time','nc_int',did_time);
    netcdf.putAtt(ncid,id_time,'standard_name','time');
    netcdf.putAtt(ncid,id_time,'units','seconds since 2001-01-01');
    netcdf.putAtt(ncid,id_time,'calendar','standard');
    netcdf.putAtt(ncid,id_time,'comment','timestamp given at the end of the sample');
    netcdf.putAtt(ncid,id_time,'ancillary_variables','sample_tms')
    if totsamp_flag
        netcdf.putAtt(ncid,id_time, 'note', 'dublicate time stamps found in the lv0-file, the first occurrence of the dublicate time is removed')
    end
    
end


function id_lat = lat(ncid)
    id_lat = netcdf.defVar(ncid,'lat','nc_float',[]);
    netcdf.putAtt(ncid,id_lat,'long_name','latitude');
    netcdf.putAtt(ncid,id_lat,'standard_name','latitude');
    netcdf.putAtt(ncid,id_lat,'units','degrees_north');
    netcdf.putAtt(ncid,id_lat,'_FillValue',NaN('single'));
end


function id_lon = lon(ncid)
    id_lon = netcdf.defVar(ncid,'lon','nc_float',[]);
    netcdf.putAtt(ncid,id_lon,'long_name','longitude');
    netcdf.putAtt(ncid,id_lon,'standard_name','longitude');
    netcdf.putAtt(ncid,id_lon,'units','degrees_east');
    netcdf.putAtt(ncid,id_lon,'_FillValue',NaN('single'));
end

function id_sampleTms = sampleTms(ncid, did_time)
    id_sampleTms = netcdf.defVar(ncid,'sample_tms','nc_int',did_time);
    netcdf.putAtt(ncid,id_sampleTms,'long_name','milliseconds of sample');
    netcdf.putAtt(ncid,id_sampleTms,'units','millisec');
    netcdf.putAtt(ncid,id_sampleTms,'_FillValue',int32(-999));
    netcdf.putAtt(ncid,id_sampleTms,'comment','add to time to get the precise sample time');  

end

function id_no_seq = chirp_sequence(fileid, did_no_seq)
    id_no_seq = netcdf.defVar(fileid,'chirp_sequence','nc_int',did_no_seq);
    netcdf.putAtt(fileid,id_no_seq,'long_name','chirp sequences');
    netcdf.putAtt(fileid,id_no_seq,'units','count');
    netcdf.putAtt(fileid,id_no_seq,'_FillValue',int32(-999));
end


function id_swv = radar_software(fileid)
    id_swv = netcdf.defVar(fileid,'radar_software','nc_float',[]);
    netcdf.putAtt(fileid,id_swv,'long_name','radar software version');
    netcdf.putAtt(fileid,id_swv,'_FillValue',NaN('single'));
end

function id_DoppMax = DoppMax(fileid, did_no_seq)
    id_DoppMax = netcdf.defVar(fileid,'nyquist_velocity','nc_float',did_no_seq);
    netcdf.putAtt(fileid,id_DoppMax,'long_name',['Nyquist velocity']);
    netcdf.putAtt(fileid,id_DoppMax,'units','m s-1');
    netcdf.putAtt(fileid,id_DoppMax,'_FillValue',NaN('single'));
    netcdf.putAtt(fileid,id_DoppMax,'comment','max. unambigious Doppler velocity; due to dealiasing values outside of the +-Nyquist velocity range are possible');
end

function id_SeqIntTime = SeqIntTime(fileid, did_no_seq)
    id_SeqIntTime = netcdf.defVar(fileid,'chirpseq_int_time','nc_float',did_no_seq);
    netcdf.putAtt(fileid,id_SeqIntTime,'long_name','integration time of each chirp sequence');
    netcdf.putAtt(fileid,id_SeqIntTime,'units','sec');
    netcdf.putAtt(fileid,id_SeqIntTime,'_FillValue',NaN('single'));
end

function id_range_offsets = range_offsets(fileid,did_no_seq)
    id_range_offsets = netcdf.defVar(fileid,'chirpseq_startix','nc_int',did_no_seq); % used to be range_offsets
    netcdf.putAtt(fileid,id_range_offsets,'long_name','chirp sequence start index');
    netcdf.putAtt(fileid,id_range_offsets,'_FillValue',int32(-999));
    netcdf.putAtt(fileid,id_range_offsets,'comment','index for the height layer where each chirp sequence starts; zero-based indexing');
end


function id_SLv = SLv(fileid, did_height, did_time)
    id_SLv = netcdf.defVar(fileid,'sensitivity','nc_float',[did_height,did_time]);
    netcdf.putAtt(fileid,id_SLv,'long_name','linear sensitivity limit');
    netcdf.putAtt(fileid,id_SLv,'units','dBz');
    netcdf.putAtt(fileid,id_SLv,'_FillValue',NaN('single'));
    netcdf.putAtt(fileid,id_SLv,'comment','provided by radar software, data only logged for bins where signal exceeds this threshold');
end