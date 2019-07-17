function [tempstruct, no_clean_signal, idx_0] = dealias_spectra_find_nonaliased_bin(cth_fin, cbh_fin, spec, range_offsets, vel, nAvg, moment_string, nf_string, nf, nbins, alias_flag, noise)

% this function looks for a range bin where no aliasing occurs. if there is
% a clean bin, the function will calculate the higher moments in this bin

% input:
%   see main function dealias_spectra.m
%
% output:
%   tempstruct: struct containing moments specified in 'moment'
%   no_clean_signal: logical == true when all bins are aliased
%   idx_0: index of range bin where the first clean signal was found

leave = false;
n_cloudbins = cth_fin-cbh_fin + 1;
cc = 0;
no_clean_signal = false;
idx_0 = NaN;

while leave == false && cc < n_cloudbins && no_clean_signal == false % look until non-aliased bin with signal after radar_moments.m is found
    % now start on top and find first column with signal that doesn't have
    % aliasing
    idx_0 = find(alias_flag(cbh_fin:cth_fin-cc) == 0 & ~isnan(spec(cbh_fin:cth_fin-cc,1)),1,'last') + cbh_fin - 1;
    if isempty(idx_0) % dealiasing not possible in the whole layer
        no_clean_signal = true;
        tempstruct = NaN;
        continue
    end
    
    % calculate moments for this range gate
    Nfft = sum(~isnan(spec(idx_0,:)));
    
    % get range indexes
    r_idx = dealias_spectra_get_range_index(range_offsets, idx_0);
    
    tempnoise.meannoise = noise.meannoise(idx_0);
    tempnoise.peaknoise = noise.peaknoise(idx_0);
    
    % calculate moments
    tempstruct = radar_moments(spec(idx_0,1:Nfft),vel(1:Nfft,r_idx),nAvg(r_idx),'noise',tempnoise,'moment_str',moment_string,'linear',nf_string,nf,'nbins',nbins);
    
    
    if isnan(tempstruct.vm) % no significant signal was detected by radar_moments.m
        % write moment to output in case they are not dealiased
        cc = cc + 1;
    else
        leave = true;
    end
end % while
