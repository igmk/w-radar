% exmple configuration file to run data processing
% 
% Configuration of the processing script. It includes input and output
% directories of the data, sets radar options (radar location, radar name),
% processing options (dealiaising, overwriting of output data) and others.
% Meta data for output file is given in a separate file.
%
% Processing is done by Raw2l_radar.m (main routine) 
% Example Raw2l_radar('config_joyrad','yyyymmdd')
%
% -------------------------------------------------------------------------
% Settings configuration file
%--------------------------------------------------------------------------

% define inpput data path - structer in selected folder should be 
%                           /yyyy/mm/dd/file-type
config.datapath = '/path/to/data/l0';

% select input-file-type
%   needs to be a case sensitive match to the file name
config.filetype = '*lv0';

% define output path - folder structure will be: /yyyy/mm/dd/*.nc
config.outputpath = '/path/to/data/l1';

% Give type of output: use one of the provided options or add your own.
% Supported:
%   ucoldefault
%   ucoldefaultcompact
%   ucoldefaultall
%   geoms
%   eurec4a
% 
% To add your own custom type, add:
%    - function outname_<config.outputtype> to determine the output file name
%      convention (see findoutfilename.m for examples), 
%    - function write_<config.outputtype> to call your write functions (see
%      savedata.m for examples
%    - your custom function(s) that write(s) the data into a file
%
% This option replaces the previous compact_flag and compact_flag_geoms.
config.outputtype = 'ucoldefault';

% All user provided metadata has been moved from this config file to another
% file - give filename here. Does not need to match with config.outputtype.
config.metadatafile = 'outputmeta_ucoldefault';


% Tag added to the file name before the .nc
config.filetag = 'v2.0_202007' ;

% Instrument nickname:
% - apears in the output file name
% - name has to match the possible preprocessing and postprocessing matlab-
%   function file names
config.nickradar = 'joyrad94';

% Station nick name (used for file naming):
% - naming convention to use: 3 letters refering to the station measured
%   example: JOYCE = joy
config.nickstation = 'joy';   

% Debuging option
% 0 - no debugging information given, code does not crash
% 1 - makes check plot in dealising, enables debugging of some functions
config.debuging = 0; 

%Overwrite existing data
% 0 - if output file(s) already exist, no data processing is done
% 1 - overwrite (existing) output file(s) 
config.overwrite = 1; 

%Dealias:
% true - dealiasing is applied
% fales - dealiasing not appleied
config.dealias = true;

%Speckle filter (optional):
% 0 - no speckle filtering applied
% 1 - simple speckle removal in time-height space
config.speckle = 1;

% moments_cal:
% Note! at the moment (July 2019), only use 2!
% 2 - means spectral moments will be calculated (runs momentslv0-function)
% 3 - means moments are taken from files lv1 if available, (runs momentslv1-function, in 201907 not working)
config.moments2calculate = 2;

% Flag to ignore LWP from the output file - optional. 
% config.ignoreLWP = 1;
