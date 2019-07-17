function [alias_flag, noise] = dealias_spectra_check_aliasing(ss, spec, vel, nAvg, range_offsets, varargin)

% varargin can contain noise

noiseflag = false;

% check and flag aliasing
alias_flag = zeros(ss(1),1);

if isempty(varargin)
    noise.peaknoise = NaN(ss(1),1);
    noise.meannoise = NaN(ss(1),1);
else
    noiseflag = true;
    noise = varargin{1};
end
    
    

for i = 1:ss(1)

    r_idx = dealias_spectra_get_range_index(range_offsets, i);
    
    if isnan(spec(i,1)) % no signal was recorded        
        continue
    end
    
    Nfft = sum(~isnan(vel(:,r_idx)));
    
    if noiseflag == false
        % get peak noise of spectrum
        tempnoise = hildebrand_sekon(spec(i,1:Nfft),nAvg(r_idx),'mean');
        noise.meannoise(i) = tempnoise.meannoise;
        noise.peaknoise(i) = tempnoise.peaknoise;
    end
    
    % check if aliasing occured by checking if more than 'frac' percent of the bins exceeded
    % mean noise leveöl at one of the spectra
    frac = 5;    
    frac = ceil(Nfft/100*frac);
    
    n_start = sum(spec(i,1:frac+1) > noise.meannoise(i));
    n_end = sum(spec(i,Nfft-frac:Nfft) > noise.meannoise(i));
    
    if n_start >= frac || n_end >= frac % then aliasing detected 
        alias_flag(i) = 1;
    end    
        
end % i