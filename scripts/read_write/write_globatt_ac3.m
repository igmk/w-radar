function write_globatt_ac3(ncid, config, data, n_file)
% function to write global attributes to the 3 output files following ac3
% data structure
% RG 16.9.2022


glob = netcdf.getConstant('NC_GLOBAL');


netcdf.putAtt(ncid,glob,'title'             ,config.title{n_file});
netcdf.putAtt(ncid,glob,'institution'       ,config.institution);
netcdf.putAtt(ncid,glob,'source'            ,config.source     );

netcdf.putAtt(ncid,glob,'comment'           ,config.comment{n_file} );

netcdf.putAtt(ncid,glob,'chirp_program'     ,[num2str(data.progno, '%02d') ' ' data.progname ]);
if n_file~=3 % housekeeping file does not follow conventions
    netcdf.putAtt(ncid,glob,'conventions'       ,config.conventions);
end
if n_file~=1 % omitting the history attribute from the moments file
    netcdf.putAtt(ncid,glob,'history'           ,config.history    );
end
netcdf.putAtt(ncid,glob,'processing_script' ,config.processing_script);
netcdf.putAtt(ncid,glob,'processing_date'   ,datestr(now));

netcdf.putAtt(ncid,glob,'featureType'       ,config.featureType );
netcdf.putAtt(ncid,glob,'references'        ,config.references );
netcdf.putAtt(ncid,glob,'license'           ,config.license );


netcdf.putAtt(ncid,glob,'PI'                ,config.pi);
netcdf.putAtt(ncid,glob,'author'            ,config.author);
netcdf.putAtt(ncid,glob,'project'           ,config.project);


% optional notes, if necessary
if data.AntiAlias == 0 % if moments not calculated by our script, add note
    netcdf.putAtt(ncid,glob,'note_dealiasing', 'No dealiasing applied');
elseif data.AntiAlias == 1 % if moments not calculated by our script, add note
    netcdf.putAtt(ncid,glob,'note_dealiasing', 'Dealising by radar software');
end

if data.CompEna > 0 % if spectra compressed in radar software, add note
    netcdf.putAtt(ncid,glob,'note_compression', 'Spectral compression enabled in RPG software');
end
