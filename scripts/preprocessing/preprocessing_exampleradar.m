function [data, config] = preprocessing_exampleradar(data, config, lv0filetype)

% this moved to read_netcdf_v1.m, fits better there
% in RPG software version 1 reflectivity was given in dB, converted to 
% linear to be compatible with further processing
% if lv0filetype == 1 % convert into linear regime
%     data.spec = 10.^(data.spec./10);    
%     data.Ze = 10.^(data.Ze/10);
% end


