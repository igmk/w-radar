function [spec_out, vel_out, status_flag, moments] =...
    dealias_spectra_from_idxA_to_idxB(idxA,...  %  start index
    idxB,...    %   end index
    range_offsets,... % index where chirp sequences start
    vel,... % velocity array containing Doppler velocity bins for each chirp sequence
    delv,... %  Doppler velocity resoliton
    spec,... % spectra full column
    vn,... % nyquist velocity for each chirp sequence
    moments,... % struct containing radar moments in the entire column
    moment_string,...  %string indicating highest moment to be calculated
    nAvg,... % number of spectral averages
    nf_string,...% indicating if mean or peak noise factor is taken for radar_moments
    nf,... % noise factor
    nbins,... % number of consecutive bins in for peak detection in radar_moments
    status_flag,...% status flag of full column
    dr,... % range resolution
    vm_prev_col,...% mean doppler velocity of previous column
    peaknoise) % peak noise of raw spectra

% output:
%   spec_out: dealiased spectra of layer
%   vel_out: velocity array for each spectrum
%   status_flag = true, if problem occured
%   moments: radar moments calculated from dealiased spectra


% ########### check into which direction the dealiasing should take place
if idxA < idxB % from bottom to top
    inc = 1;
else
    inc = -1;
end

% ######### preallocate variables
ss = size(spec);
idx = idxA:inc:idxB;
n_levels = numel(idx); % from top to bottom

spec_out = spec(idx,:);
vel_out = NaN(n_levels,ss(2));

tempflag = false(size(ss(1))); % indicates if spectra were determined correctly

Nfft = sum(~isnan(vel));

cc = 0;
for ii = idxA:inc:idxB % 
    cc = cc + 1;
    

    % ############### get range indexes
    r_idx = dealias_spectra_get_range_index(range_offsets, ii);
    
    if any( isnan(spec(ii,1:Nfft(r_idx))) ) || sum(spec(ii,1:Nfft(r_idx))) < 10^-20 % then no signal is available
        continue
    end
    
%     if ii == 227
%         ii
%     end
    % ############## get spec and vel chains
    
    vel_chain = (vel(1,r_idx) - 2*Nfft(r_idx)*delv(r_idx)) : delv(r_idx) : (vel(Nfft(r_idx),r_idx)+2*Nfft(r_idx)*delv(r_idx) + delv(r_idx));    
    
    % check if chrip boundary is approached
    [~,idx] = min(abs(ii-range_offsets));
    next_chirp = range_offsets(idx);
    
    % check if boundaries are exceeded
    if ii - inc > numel(moments.vm) || ii -inc < 1
        vm_guess = NaN;
    else
        vm_guess = moments.vm(ii-inc);
    end
    
    % check here if vm_guess has signal
    if isnan(vm_guess)
        vm_guess = dealias_spectra_vm_guess_qual_check(moments.vm, vm_prev_col, ii, inc, dr(r_idx));
    end
    
    % check if upper or lower boundary is reached
    [spec_chain, status_flag(ii,1:4)] = dealias_spectra_concetenate_spectra(vm_guess, spec(:,1:Nfft(r_idx)), vn(r_idx), ii, next_chirp, Nfft(r_idx));
    
    
    %################ get final spectrum    
    [spec_out(cc,1:Nfft(r_idx)), vel_out(cc,1:Nfft(r_idx)), status_flag(ii,1:4)] = dealias_spectra_determine_final_spectrum(vm_guess, spec_chain, vel_chain, Nfft(r_idx));
    
    
    %############## quality check final spectrum
    alias_flag = dealias_spectra_quality_check_final_spectrum(spec_out(cc,1:Nfft(r_idx)), peaknoise(ii));
    
  
    % ############ if alias_flag == 1, then try centering again with
    % modified vm_guess and get vm
    if alias_flag == 1
        % use adjacent velocities to get a new vm_guess, make sure that
        % indexes do not exceed boundaries
        a = ii - 2; if a == -1, a = 1; elseif a == 0, a = 2; end
        b = ii + 2; if b == ss(1) + 2, b = ss(1); elseif b == ss(1) + 1, b = ss(1) - 1;  end
        
        vm_guess_new = dealias_spectra_vm_guess_modify(vm_prev_col(a:b), vn(r_idx), vm_guess);

        if ne(vm_guess_new, vm_guess) && ~isnan(vm_guess_new) % then determine spectrum again
            [spec_out(cc,1:Nfft(r_idx)), vel_out(cc,1:Nfft(r_idx)), status_flag(ii,1:4)] = dealias_spectra_determine_final_spectrum(vm_guess_new, spec_chain, vel_chain, Nfft(r_idx));
            
            %############## quality check final spectrum again
            alias_flag = dealias_spectra_quality_check_final_spectrum(spec_out(cc,1:Nfft(r_idx)), peaknoise(ii));
        end
               
    end
   
  
    % ############## calculate moments
    tempstruct = radar_moments(spec_out(cc,1:Nfft(r_idx)),vel_out(cc,1:Nfft(r_idx)),nAvg(r_idx),'moment_str',moment_string,'linear',nf_string,nf,'nbins',nbins);
    moments = dealias_spectra_write_tempmoments_to_finalmoments(moments, tempstruct, ii, moment_string);

    
    if alias_flag == 1 % then centering did not work properly all follwing bins might be affected
        tempflag(ii:inc:idxB) = true;
    end
        
    
    
    
end % ii


if idxA > idxB % flip order of output
    spec_out = flip(spec_out,1);
    vel_out = flip(vel_out,1);
end

% update status_flag
status_flag(tempflag == true,2) = '1';
