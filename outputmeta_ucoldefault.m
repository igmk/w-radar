% Contains user defined meta data for the netcdf file(s)
% produced. The name of this file needs to correspond
% to the filename given in the config_***.m file
%
% Meta data required for the University of Cologne   
% default ouput files.


% height above mean sea level [m] of the instrument platform
%   stored in the metadata 
config.MSL = 11;

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
% Name of the processing script:
config.processing_script = 'https://github.com/igmk/w-radar/releases/tag/v2.0_202007';

