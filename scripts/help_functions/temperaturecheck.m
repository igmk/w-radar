function data = temperaturecheck(data)
% Sometimes in the temperature variables found 0 values. This function 
% replacing these with NaNs since a temperature of 0 K is obviously wrong
% RG 20.10.2022


data.T_env( data.T_env == 0 ) = NaN;
data.T_pc( data.T_pc == 0 ) = NaN;
data.T_rec( data.T_rec == 0 ) = NaN;
data.T_trans( data.T_trans == 0 ) = NaN;