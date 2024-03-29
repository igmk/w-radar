function [] = savedata(data, config)
% Updated to be flexible for different output filenaming conventions.
% RG 13.4.2022

% finding output filename moved to it's own function findoutfilename - RG


% if output file already exists (could happen when overwrite requested, or 
% when only some of the requested output files exist), delete the existing
% files


fn = fieldnames(config.outfiles);
for k=1:numel(fn) % loop over output files
    
    if exist(config.outfiles.(fn{k}), 'file')
        delete(config.outfiles.(fn{k}));
    end

end

% create directories, in case missing (currently multiple output
% directories not supported)

[pathstr,~ ,~] = fileparts(config.outfiles.(fn{1}));
if ~exist(pathstr,'dir')
    mkdir(pathstr);
end


% call write function(s) for output type
outfilefunc = str2func(['write_' config.outputtype]);
outfilefunc(data, config);



end % function




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% local functions for calling write functions
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


function write_ucoldefaultall(data, config)
    % write output file for U. Cologne default option (full file)
    % RG 13.4.2022
    
    disp(['Writing file ' config.outfiles.file1])
    write_joyrad94_data_2_nc(data,config.outfiles.file1, config);  
    
end % function


function write_ucoldefaultcompact(data, config)
    % write output file for U. Cologne default option (compact file)
    % RG 13.4.2022
    
    disp(['Writing file ' config.outfiles.file1])
    write_joyrad94_data_2_nc_compact(data,config.outfiles.file1, config);
    
end % function


function write_ucoldefault(data, config)
    % write output file for U. Cologne default option (compact and full file)
    % RG 13.4.2022
    
    disp(['Writing file ' config.outfiles.file1])
    write_joyrad94_data_2_nc(data,config.outfiles.file1, config); 
    
    disp(['Writing file ' config.outfiles.file2])
    write_joyrad94_data_2_nc_compact(data,config.outfiles.file2, config);
    
end % function


function write_geoms(data, config)
    % write output file for geoms-type output
    % RG 13.4.2022
    
    disp(['Writing file ' config.outfiles.file1])
    write_joyrad94_data_2_nc_geoms(data,config.outfiles.file1, config);
    
end % function


function write_eurec4a(data, config)
    % write output file for EURECa campaign: physical and technical
    % parameters file
    % RG 13.4.2022
    
    % calculating Doppler velocity arrays for output
    data.velocitymatrix = calculate_doppler_arrays(data);
    
    disp(['Writing file ' config.outfiles.file1])
    write_data_2_nc_physparam(data, config.outfiles.file1, config);
        
    disp(['Writing file ' config.outfiles.file2])
    write_data_2_nc_technical(data,config.outfiles.file2, config);
     
    
end % function



function write_ac3(data, config)
    % write output files for ac3 data set
    % RG 1.6.2022
    
    disp(['Writing file ' config.outfiles.file1])
    write_data_2_nc_moments(data, config.outfiles.file1, config);
        
    disp(['Writing file ' config.outfiles.file2])
    write_data_2_nc_spectra(data, config.outfiles.file2, config);
     
    disp(['Writing file ' config.outfiles.file3])
    write_data_2_nc_housekeep(data, config.outfiles.file3, config);
    
end % function
