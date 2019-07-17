function alias_flag = dealias_spectra_quality_check_final_spectrum(spec, peaknoise)

% check if aliasing occured
% a certain fraction of all values that exceed the peak noise level must be
% located within +-v_n/4, otherwise the spectrum was not aligned correctly

% input:
%   spec: linear Doppler spectrum (1 x Nfft)
%   peaknoise: linear peaknoise determined before dealiasing
%
% output:
%   alias_flag == 1 if aliasing was not performed correctly

alias_flag = 0;

ss = size(spec);


% sort spectrum
[spec_sorted, idx] = sort(spec);

% get all values above peaknoise
idx_signal = spec_sorted > peaknoise;

% select only the highest 50%
n_signal = sum(idx_signal);
idx_signal(1:end-ceil(n_signal/2)) = false;

% check how many of these points are within v(1)/v(end) -+v_n/4
idx_inside = idx(idx_signal) < ss(2)/4 | idx(idx_signal) > 3/4*ss(2);


% calculate fraction of how mmany percents fall into v(1)/v(end) -+v_n/4
frac = sum(idx_inside)/sum(idx_signal);

if frac > 0.8 % then more than 80 % of the largest 50 % are located at the edge
    alias_flag = 1;
end



