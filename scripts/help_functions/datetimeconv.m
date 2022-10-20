
function timeout =  datetimeconv(yy,mm,dd,varargin) % datetimeconv(yy,mm,dd,hh,mmm,ss)
% function to convert date and time to the time format used in the program
% RG 9.7.2019

 if isempty(varargin)
     hh = 0;
     mmm = 0;
     ss = 0;
 else
     hh = varargin{1};
     mmm = varargin{2};
     ss = varargin{3};
 end

    timeout  = datenum([yy,mm,dd,hh,mmm,ss])*86400 - datenum([2001 1 1])*86400;
    
end % function

