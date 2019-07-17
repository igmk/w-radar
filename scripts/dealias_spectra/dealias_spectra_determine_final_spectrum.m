function [spec_out, vel_out, status_flag] = dealias_spectra_determine_final_spectrum(vm_guess, spec_chain, vel_chain, Nfft)

% function determines the final spectrum by
%   1) using vm_guess as initial guess for location of the true spectrum's
%   peak
%   2) find the maximum that is closes to vm_guess
%   3) center the spectrum around the maximum so that a noise floor is
%   present at the left and right edge of the spectrum
%
% input:
%   vm_guess: intial guess velocity [m/s]
%   spec_chain: 5 concenated spectra, spectrum #3 is the considered bin
%   vel_chain: velocity array corresponding to spec_chain
%   Nfft: number of bins per spectra
%
% output:
%   spec_out: final spectrum (1 x Nfft)
%   vel_out: final velocity array (1 x Nfft)
%   status_flag == true if a clean determination could not be performed
%       e.g.:
%           - vm_guess == NaN
%           - the peak of the final spectrum was close to the boundaries of
%               vel_chain


status_flag = '0000';

% ############### find maximum of spectrum

if isnan(vm_guess) % start at center of 3rd spectrum
    idx_in = single(floor(5/2*Nfft));
    status_flag(4) = '1';
else
    % find index nearest value in vel_chain to vm_guess
    [~,idx_in] = min(abs(vm_guess-vel_chain));
end

% look in the area of idx_in +-Nfft/4 for the maximum
% check if boundary is exceeded
    
if idx_in-Nfft/4 < 1
    [~,idx_max] = max(spec_chain(1:Nfft));
    status_flag(3) = '1';
elseif idx_in+Nfft/4 > numel(spec_chain)
    [~,idx_max] = max(spec_chain(end-Nfft:end));
    idx_max = idx_max + numel(spec_chain) - Nfft - 1;
    status_flag(3) = '1';
else
    [~,idx_max] = max(spec_chain(idx_in-Nfft/4:idx_in+Nfft/4));
    idx_max = idx_max + idx_in-Nfft/4 - 1;
end
    
    


% ################ if the spectrum is broad then the new spectrum can still
% contain aliased contributions. to minimize that influence
% cetner the spectrum so that it has the lowest value when adding
% the signals of the first and last entry, respectively.

% check first that boundaries are not crossed
if idx_max - Nfft/4 - Nfft/2 >= 1 && idx_max + Nfft/4 + Nfft/2 - 1 <= numel(spec_chain)
    
    % then neither lower nor upper index (i.e. 1 or numel(spec_chain))
    % will be exceeded by the following procedure.
    shift_factor = 4;
    
    idx_shift = -Nfft/shift_factor:1:Nfft/shift_factor;
    idx_edges = [(idx_max-3*Nfft/shift_factor:idx_max-Nfft/shift_factor)', (idx_max+Nfft/shift_factor-1:idx_max+3*Nfft/shift_factor-1)'];
    edgesum = sum(spec_chain(idx_edges),2);
    
    % get minimum of edgesum
    [~,idx_min] = min(edgesum);
    idx_max = idx_shift(idx_min) + idx_max;
    
end



% ################# create new spectrum

if idx_max-Nfft/2 < 1
    spec_out = spec_chain(1:Nfft);
    vel_out = vel_chain(1:Nfft);
elseif idx_max+Nfft/2-1 > numel(spec_chain)
    spec_out = spec_chain(end-Nfft+1:end);
    vel_out = vel_chain(end-Nfft+1:end);
else
    spec_out = spec_chain(idx_max-Nfft/2:idx_max+Nfft/2-1);
    vel_out = vel_chain(idx_max-Nfft/2:idx_max+Nfft/2-1);
end


