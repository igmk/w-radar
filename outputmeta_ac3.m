% Contains user defined meta data for the netcdf file(s)
% produced. The name of this file needs to correspond
% to the filename given in the config_***.m file
%
% Meta data required for the AC3 ouput files.


% height above mean sea level [m] of the instrument platform
%   stored in the metadata 
config.MSL = 0;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Standard global attributes following the CF convention (descriptions from CF documentation)
config.conventions = 'CF-1.10';

% "A succinct description of what is in the dataset."
config.title = {'Radar moments', ... % string for title attribute of each file
                'Doppler spectra', ...
                'Housekeeping data'};

% "Provides an audit trail for modifications to the original data. Well-behaved generic netCDF filters will automatically append their name and the parameters with which they were invoked to the global history attribute of an input netCDF file. We recommend that each line begin with a timestamp indicating the date and time of day that the program was executed."
config.history  = [datestr(now) ': Lv1 file generated from binary lv0 file using processing_script'];

% "Specifies where the original data was produced."
config.institution = '';

% "The method of production of the original data."
config.source = '';%'RPG-FMCW-94-SP';

% "Miscellaneous information about the data or methods used to produce it."
% separate comments for each of the three files, 
config.comment = {'This file contains radar moments and the most important parameters needed to analyze the data. Doppler spectra and full housekeeping data are stored in separate files.', ... % moments file
                  'The Doppler spectra has gone through a dealiasing algorithm, for details see references. This file only contains the radar Doppler spectra and associated parameters. Further radar variables and full housekeeping data are stored in separate files.', ... % spectra file
                  'This file only contains housekeeping data, 89 GHz brightness temperature and derived products, and weather station measurements. Main radar variables are stored in separate files.'};    % housekeeping file

% "Published or web-based references that describe the data or methods used to produce it."
config.references  = '';


% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Other meta data

% Name of the processing script
config.processing_script = '';

config.featureType = 'timeSeriesProfile';

config.pi          = '';       
config.author      = '';       
config.project     = '';

config.license     = '';
