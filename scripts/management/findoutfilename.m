function [config] = findoutfilename(config, filein)
% Function to find output file name. Moved here from savedata.
% RG 27.6.2019
% 
% Updated to be flexible for different output filenaming conventions.
% RG 13.4.2022

[~, filename, ~] = fileparts(filein); 
config.filename = filename;

%Remove from filename duplicated info
filename = strrep(config.filename,config.nickradar, '');
filename = strrep(filename,config.nickstation, '');

if isempty( strfind(config.filetag, '_') ) && ~isempty(config.filetag)

    config.filetag = [ '_' config.filetag];
end


% Get output filenames
outfilefunc = str2func(['outname_' config.outputtype]);
outfile = outfilefunc(config, filename);


% loop over all outfiles (number varies depending on set-up)
fn = fieldnames(outfile);
for k=1:numel(fn)
    
    % remove possible doublicate underscore
    while strfind(outfile.(fn{k}), '__')
        outfile.(fn{k}) = strrep(outfile.(fn{k}),'__', '_');
    end

end

config.outfiles = outfile;


end % function



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% local functions for output file names
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function out = outname_ucoldefaultall(config, filename)
    % output filename for U. Cologne default option (full file)
    % RG 13.4.2022

    out.file1 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
    
end % function


function out = outname_ucoldefaultcompact(config, filename)
    % output filename for U. Cologne default option (compact file)
    % RG 13.4.2022

    out.file1 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s_compact%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
    
end % function


function out = outname_ucoldefault(config, filename)
    % output filename for U. Cologne default option (compact and full file)
    % RG 13.4.2022

    out.file1 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
    out.file2 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s_compact%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
    
end % function

function out = outname_geoms(config, filename)
    % filename for geoms-type output
    % RG 13.4.2022

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Modified by A.Pschera on 2021-11-17 - added different naming for
    %mirac:
    if config.nickradar([1:5]) == 'mirac'
        out.file1 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s_%s_%s_%s%sT%s*.nc','groundbased','radar_profiler',config.nickradar,config.file_association,config.nickstation,'20',filename([6:11]),filename([12:17])));
    else
        out.file1 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s_%s_%s_%s%sT%s*.nc','groundbased','radar_profiler',config.nickradar,config.file_association,config.nickstation,'20',filename([1:6]),filename([8:13])));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end % function

function out = outname_eurec4a(config, filename)
    % output filename for W-band radar on MS Merriam during the EUREC4A campaign
    % RG 13.4.2022

    out.file1 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
    out.file2 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s_technical%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
    
end % function




function out = outname_ac3(config, filename)
    % output filenames used for ac3 data set
    % RG 1.6.2022
    
    stripped_filestr = strrep(filename,'_ZEN', '');
    stripped_filestr = strrep(stripped_filestr,'ZEN', '');

    out.file1 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s_%s.nc', config.nickradar, config.nickstation, config.filetag, stripped_filestr));
    out.file2 = fullfile(config.outputpath_tree, sprintf('%s_%s_spectra_%s_%s.nc', config.nickradar, config.nickstation, config.filetag, stripped_filestr));
    out.file3 = fullfile(config.outputpath_tree, sprintf('%s_%s_housekeep_%s_%s.nc', config.nickradar, config.nickstation, config.filetag, stripped_filestr));
    
end % function

