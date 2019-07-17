
This Matlab program reads data recorded with the FMCW-radar of RPG. RPG 
stores the data in binary format. The program converts it into netcdf4 
format. Additionally some reprocessing of the data can be applied including 
dealasing the radar spectra and the calculation of higher moments. 

The program was first created by Nils Küchler in 2018. The code was 
restructured and partly updated in 2019 by Juan Antonio Bravo Aranda and 
Rosa Gierens.


## Content of this documentation: ##
0. General stucture
1. How to call the program
2. Setting data processing options (config file)
3. Input and output files

    3.1. Supported input files
    
    3.2. output files: naming convention and types
    
4. Short description of the data processing

    4.1 Steps of the data processing

    4.2 Dealiasing


####################### DESCRIPTION OF THE PROCESSING #####################


# 0. General stucture #

All code for this program is contained in direcotry scripts, and its 
subdirectories. It is assuemd that "scripts" is located in the same 
directory where you execute the program. The config file should also be in
the same directory, otherwise give full path of the config file instead of
config_radarname (see below).


# 1. How to call the program #

Call the scripts with:
- Raw2l1_radar('config_radarname')
    to process today's and yesterdays data
- Raw2l1_radar('config_radarname', 'yyyymmdd')
    to process the given date
- Raw2l1_radar('config_radarname', 'yyyymmd1', 'yyyymmd2') 
    to process all data from the first date to the second date



# 2. Setting data processing options (config file) #

The settings required for running the program have to be set in the 
configuration file (config_radarname, see above for calling the program). 
This file includes:
- input and output directories of the data
- sets radar options (radar location, radar name, contact person)
- processing options:
    - dealiasing
    - overwriting of output file(s)
    - debugging
    - definition of output file type
    - whether moments are calculated from spectra or copied from RPG lv1 
      file (this option is not available at the moment 10.7.2019)

See config_example.m for details.


# 3. Input and output files #

## 3.1 Supported input files ##

The program can handle 4 different file formats for lv0 files:
    i) binary files of type 1 created with RPG radar software version 1.
    ii) netcdf files created from i) where no additional processing 
        has been applied. 
    iii) binary files of type 2 created with RPG radar software version 2. 
    iv) binary files of type 3 created with RPG radar software version 3-5.
This program automatically identifies the file type and adjusts reprocessing
and creates a unified netcdf file for any of the file types i) to iv). 
Input data directory defined in config file (config.datapath). Data is 
assumed to be located in config.datapath/yyyy/mm/dd/
Note that you need to have read access to the directory.
The option of reading RPG lv1 binary files and converting them to netcdf is
not working at the moment (10.7.2019).

## 3.2 Output files: Naming convention and types ##

Two types of putput files are available. 
    i) general file: includes all metadata information, all flags, 
       all spectra, all moments
    ii) compact file: includes only moments (no spectra), some metadata 
Set in the config file if i), ii), or both should be created.

Naming convention for output files:
i) radarname_station_yyyymmddHHMMSS_program_scan.nc,
    example: joyrad94_nya_20190710000000_P01_ZEN.nc
ii) radarname_station_yyyymmddHHMMSS_program_scan_compact.nc,
    example: joyrad94_nya_20190710000000_P01_ZEN_compact.nc
where radarname and station are specified in config file, program refers to 
the program number (chirp table) defined in the RPG radar software, and 
scan is the scanning strategy also set in the RPW radar software. 'program' 
and 'scan' are copied from the input file name.

Output data directory defined in config file (config.outputpath). Data will
be written in config.outputpath/yyyy/mm/dd/*.nc If directories are not 
existing they will be created by the program. Note that you need to have 
write access to the target directory.


# 4. Short description of the data processing (Raw2l1_radar) #

Raw2l1_radar includes:
- Configures the time period to be processed
- Load config information
- Run code for each day selected: momentslv0
    - In momentslv0.m first search for input files for the day considered.
    - Per selected day, ech file is processed individually. Details below.

## 4.1. Steps of the data processing ##

This is the processing applied for each file in momentslv0.m.

i) Check if output file already exist. 
    If config.overwrite is set to 0, and the desired output files already 
    exist, no processing of that file will be done. Program continues with 
    next file.

ii) Reading input data. 
    First the lv0 file type is determined, and the corresponding reading 
    function is assigned, for example read_lv0_v3.m (function whichReader.m).
    If a new file type is introduced by RPG, a new reading function needs 
    to be created. The input file is read into variable "data" (variable 
    type structure).

iii) Additional variables added into "data" and missing values set to NaN.
    Function setting_data.m. If you want to add any variables into "data", 
    do it here.

iv) Radar specific settings to be applied before dealising. 
    The function preprocessing_radarname.m is called if such a function
    exist. Such a function is only necessary if you want to make any radar 
    specific corrections or add missing meta data (relevant for early RPG 
    software versions). Radar specific corrections are for example removing 
    known artifacts from the spectra. Note that you can define corrections 
    to only be applied for specific time periods. Each radar has its own 
    function, and any adjustements are radar specific. The file name has to 
    match the config.nickradar given in the config file. For a simple 
    example of the preprocessing function see preprocessing_exampleradar.m.

v) Dealising and calculating moments 
    If set in config file, and supported by the program, the dealising 
    routine is called (function dealising.m). If dealising is not applied, 
    moments are calculated from the non-dealiased spectra (function 
    moments_retrieval.m). Details of the dealiasing in Section 4.2.

vi) Radar specific settings to be applied after dealising. 
    The function postprocessing_radarname.m is called if such a function
    exist. Such a function is only necessary if you want to make any radar 
    specific corrections, for example add known reflectivity offsets due
    calibration issues. Note that you can define corrections to only be 
    applied for specific time periods. Each radar has its own function, and 
    any adjustements are radar specific. The file name has to match the 
    config.nickradar given in the config file. For a simple example of the 
    postprocessing function see postprocessing_exampleradar.m.

vii) Write output data
    Output file(s) are created and data is written to file(s). Function
    savedata.m calls function write_joyrad94_data_2_nc.m for general file 
    and function write_joyrad94_data_2_nc_compact.m for compact file. To 
    add variables in output, edit the write***.m function(s).


######################## END OF THE DESCRIPTION ###########################
   



if called: calculates higher spectral moments calling "radar_moments.m" 
and/or dealiases radar Doppler spectra calling "dealias_spectra.m"

## dealias_spectra.m
if terminal velocity of particles exceed the nyquist limit, then, in 
FMCW radar, their signal is folded into the upper and lower range gate. by 
concetenating adjacent spectra, the original spectrum can be recovered. 
the goodness of the dealiasing procedure is stored for each bin (see further 
description of variables in netcdf files, or write_joyrad94_data_2_nc)

## radar_moments.m:
calculates higher spectral moments by

i) determining the mean and peak noise level using Hildebran-Sekon procedure

ii) determining signficant signal as blocks with at least three consecutive bins above peak noise level. if a block is found it is extended by including all adjacent bins until the signal drops the first time below the mean noise level.

iii) the mean noise floor is subtracted from valid signal and spectral moments are calculated. note that all moments are calculated for the entire spectra and not just for the main peak
