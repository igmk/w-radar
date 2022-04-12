function [config] = findoutfilename(config, filein)
% Function to find output file name. Moved here from savedata.
% RG 27.6.2019

[~, filename, ~] = fileparts(filein); 
config.filename = filename;

%Remove from filename duplicated info
filename = strrep(config.filename,config.nickradar, '');
filename = strrep(filename,config.nickstation, '');

if isempty( strfind(config.filetag, '_') ) && ~isempty(config.filetag)

    config.filetag = [ '_' config.filetag];
end

switch config.compact_flag
    case 0
        outfile = fullfile(config.outputpath_tree, sprintf('%s_%s_%s%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
    case 1
        outfile2 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s_compact%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
    case 2
        outfile = fullfile(config.outputpath_tree, sprintf('%s_%s_%s%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
        outfile2 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s_compact%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
    case 3
        outfile = fullfile(config.outputpath_tree, sprintf('%s_%s_%s%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
        outfile2 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s_technical%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
end

switch config.compact_flag_geoms
    case 0
        outfile2 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s%s.nc', config.nickradar, config.nickstation, filename, config.filetag));
    case 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Modified by A.Pschera on 2021-11-17 - added different naming for
        %mirac:
        if config.nickradar([1:5]) == 'mirac'
            outfile2 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s_%s_%s_%s%sT%s*.nc','groundbased','radar_profiler',config.nickradar,config.file_association,config.nickstation,'20',filename([6:11]),filename([12:17])));
        else
            outfile2 = fullfile(config.outputpath_tree, sprintf('%s_%s_%s_%s_%s_%s%sT%s*.nc','groundbased','radar_profiler',config.nickradar,config.file_association,config.nickstation,'20',filename([1:6]),filename([8:13])));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if exist('outfile', 'var')

    while strfind(outfile, '__')
        outfile = strrep(outfile,'__', '_');
    end

    config.outfile = outfile;
end



if exist('outfile2', 'var')
    
    while strfind(outfile2, '__')
        outfile2 = strrep(outfile2,'__', '_');
    end

    config.outfile2 = outfile2;
end
