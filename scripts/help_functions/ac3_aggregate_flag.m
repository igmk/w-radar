function flag_aggregate = ac3_aggregate_flag(data)
% function to create an aggregate quality flag for output file, that
% combines several flags in a simplified way
% RG 8.12.2022

flag_aggregate = zeros(1,data.totsamp); %

% Bit0: time
data.timeshift( isnan(data.timeshift) ) = 0; % treat missing flag as no flag set
ind = data.timeshift ~= 0 ; % time stamp has been edited
flag_aggregate(ind) = flag_aggregate(ind) + 2^0;

% Bit1: dealiasing performed
ind = any(data.Aliasmask,2); % aliasing detected in column
flag_aggregate(ind) = flag_aggregate(ind) + 2^1;

% Bit2: any of the dealising flags
data.AliasStatus( isnan(data.AliasStatus) ) = 0; % treat missing flag as no flag set
ind = any(data.AliasStatus,2); % problems with dealiasing occurred
flag_aggregate(ind) = flag_aggregate(ind) + 2^2;

% Bit3:  any of the RPG flags
data.QF( isnan(data.QF) ) = 0; % treat missing flag as no flag set
ind = data.QF ~= 0 ;
flag_aggregate(ind) = flag_aggregate(ind) + 2^3;