% exmple configuration file to run data processing
% 
% Configuration of the processing script. It includes input and output
% directories of the data, sets radar options (radar location, radar name, 
% contact person), processing options (dealiaising, overwriting of output 
% data) and others.
%
% Processing is done by Raw2l_radar.m (main routine) 
% Example Raw2l_radar('config_joyrad','yyyymmdd')
%
% -------------------------------------------------------------------------
% Settings configuration file
%--------------------------------------------------------------------------

% define inpput data path - structer in selected folder should be 
%                           /yyyy/mm/dd/file-type
config.datapath = '/data/obs/site/nya/joyrad94/l0';

% select input-file-type
%   his needs to be a case sensitive match to the file name
%config.filetype = '*nc'; % (older version)
config.filetype = '*lv0'; % (newer version)

% define output path - folder structure will be: /yyyy/mm/dd/*.nc
config.outputpath = '/data/obs/site/nya/joyrad94/l1';

% tag added to the file name before the .nc
config.filetag = 'v2.0_202007' ;

% Instrument nickname:
% - apears in the output file name
% - name has to match the possible preprocessing and postprocessing matlab-
%   function file names
config.nickradar = 'joyrad94';

% station nick name (used for file naming):
% - naming convention to use: 3 letters refering to the station measured
%   example: JOYCE = joy
config.nickstation = 'joy';   

% height above mean sea level [m] of the instrument platform
%   stored in the metadata 
config.MSL = 11;

% further info for output file (metadata for the netcdf-output file) 
%   These inforamtion will be placed in the golbal attributes of the NetCDF
%   file
% Name of the PROJECT_PI
config.pi_name           = 'Surname;Name';                                 
 % AFFILIATION of the PI
config.pi_affiliation    = 'Example University (UNI);State';               
% ADDRESS of the PIs AFFILIATION
config.pi_address        = 'Example University,Institute;Street;Postcode;City;State'; 
% Mail address of the PI
config.pi_mail           = 'pi.example@example-uni.com';                   
% Name of the organisation responsible for quality controll of the data 
config.do_name           = 'Example Institution (INST);State';             
 % officiel AFFILIATION of the DO_NAMEs\
config.do_affiliation    = 'Example Institution;State';
% ADDRESS of the DOs AFFILIATION
config.do_address        = 'Example Institution;Adress'; 
% Mail address of the DO
config.do_mail           = 'do.example@example-uni.com'; 
% Name of the responible data submitter 
config.ds_name           = 'Example Institution (INST);State'; 
% officiel AFFILIATION of the DS_NAMEs\
config.ds_affiliation    = 'Example Institution;State'; 
 % ADDRESS of the DSs AFFILIATION
config.ds_address        = 'Example Institution;Adress';
% Mail address of the DS_name
config.ds_mail           = 'ds.example@example-uni.com'; 
% Optional – brief description of the file’s containing data
config.data_description  = 'daily radar measurements at station XYZ;Country'; 
% Describes field of research to which the data  belongs and the data 
%   acquisition method
config.data_discipline   = 'Atmospheric.Physics;Remote.Sensing;Radar.Profiler'; 
% Specifies the origin of the data (EXPERIMENTAL, MODAL, or both) and the 
%   spatial characteristic of the data set.Spatial dimensions are: 
%   0D = SCALAR; 1D = PROFILE; 2D for or more dimensions = FIELD. 
%   STATIONARY for fixed locations and MOVING for moving platforms
config.data_group        = 'Experimantal;Profile;Stationary'; 
% Contains the identification of the location of the reported geophysical 
%   quantities
config.data_location     = 'Research station XYZ;Country'; 
% Consist of two information, first, the instrument type used for 
%   measurements, second, the acronym of the operating institution. It may 
%   differ from the PIs or DOs affiliation
config.data_source       = strcat('Radar.Standard.Moments.Ldr_',config.nickradar,';run by UNI');
% Name of the processing script inclouding version or the repository where 
%   the script can be downloaded:
config.processing_script = 'https://github.com/igmk/w-radar/releases/tag/v2.0_202007';

% Debuging option
% 0 - no debugging information given, code does not crash
% 1 - makes check plot in dealising, enables debugging of some functions
config.debuging = 0; 

%Overwrite existing data
% 0 - if output file(s) already exist, no data processing is done
% 1 - overwrite (existing) output file(s) 
config.overwrite = 0; 

% compact-flag: 
% 0 - only general file is created (all metadata information, all flaggs, 
%     all spectra, all moments)
% 1 - only compact file is generated (only moments, some metadata 
%     inforation)
% 2 - both files are gerenrated (both files above are generated)
% 3 - all variables are included in output, but data are split in 2 files (physical variables and technical parameters)
config.compact_flag = 2;

%%%%%%%%%%%%%%
%Added on 2021-11-10 by apschera:
% compact-flag_ geoms_naming convention: 
% 0 - old university of Cologne naming convertion for the compact file is
%     used
% 1 - the compact file name follows the GEOMS naming convertion
config.compact_flag_geoms = 1;
%%%%%%%%%%%%%%

%Dealias:
% true - dealiasing is applied
% fales - dealiasing not appleied
config.dealias = true;

% moments_cal:
% Note! at the moment (July 2019), only use 2!
% 2 - means spectral moments will be calculated (runs momentslv0-function)
% 3 - means moments are taken from files lv1 if available, (runs momentslv1-function, in 201907 not working)
config.moments2calculate = 2;
