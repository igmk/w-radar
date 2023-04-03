function data = read_lv0_v3(infile)

% reads binary files from lv0 of joyrad94 and return data in a struct
% Author: Nils Küchler
% created: 9 February 2017
% modified: 9 Feburary 2017, Nils Küchler 

%%%%%%%%%%%%%%%% open file
    fid = fopen(infile, 'r', 'l');
    data.readerror = true; % flag for error in reading file
    
    if fid == -1
        disp(['>>> error opening' infile])
        return
    end
    
    %%%%%%%%%%%%%%% read header information %%%%%%%%%%%%%%%%
    
    data.filecode = int32(fread(fid,1,'int')); % lv1-file code
    data.headerlen = int32(fread(fid,1,'int')); % header length in bytes (not including headerlen)
    data.starttime = int32(fread(fid,1,'uint')); % time of first sample in file
    data.endtime = int32(fread(fid,1,'uint')); % time of last sample in file
    data.progno = int32(fread(fid,1,'int')); % program number, as definded in the host-pc software
    data.modelno = int32(fread(fid,1,'int')); % =0 singel pol., = 1 dual pol, 2 = dual pol. LDR configuration
        % modelno contains same information as DualPol, but different
        % definition - don't use to avoid confusion. RG 19.9.2022
        
    cc = 0;
    count = 1;
    while cc == 0
        ch(count) = fread(fid,1,'char*1');
        if ch(count) == 0
            data.progname = char(ch); % null terminated string of chirp program name
            cc = 1;
        else
            count = count + 1;
        end
    end
    clear ch
    
    cc = 0;
    count = 1;
    while cc == 0
        ch(count) = fread(fid,1,'char*1');
        if ch(count) == 0
            data.custname = char(ch); cc = 1; % null terminated string of custumor name
        else
            count = count + 1;
        end
    end
    clear ch cc count
    
    data.freq = single(fread(fid,1,'single')); % radar frequency [GHz]
    data.AntSep = single(fread(fid,1,'single')); % antenna separation [m]
    data.AntDia = single(fread(fid,1,'single')); % antenna diameter [m]
    data.AntG = single(fread(fid,1,'single')); % linear antenna gain
    data.HPBW = single(fread(fid,1,'single')); % half power beam width [°]
    data.C = single(fread(fid,1,'single')); % radar constant, defined by eq. (2.1.5) radar_manual_v2
    data.DualPol = int8(fread(fid,1,'char*1')); % 0 = single pol radar, 1 = dual pol radar LDR conf., 2 = dual pol radar STSR mode
    data.CompEna = int8(fread(fid,1,'char*1')); % spectral compression flag: 0 = not compressed, 1 = spectra compressed, 2 = spectra compressed and spectral polarimetric variables stored in the file
    data.AntiAlias = int8(fread(fid,1,'char*1')); % 0 = Doppler spectra are not anti-aliased, 1 = doppler spectra have been anti-aliased
    data.SampDur = single(fread(fid,1,'single')); % sample duration [sec]
    data.Lat = single(fread(fid,1,'single')); % GPS latitude
    data.Lon = single(fread(fid,1,'single')); % GPS longitude
    data.CalInt = int32(fread(fid,1,'int')); % period of automatic zero calibrations in number of samples
    data.n_levels = int32(fread(fid,1,'int')); % number of radar altitude layers
    data.T_altcount = int32(fread(fid,1,'int')); % number of temperature profile altitude layers
    data.H_altcount = int32(fread(fid,1,'int')); % number of humidity profile altitude layers
    data.no_chirp_seq = fread(fid,1,'int'); % number of chirp sequences
    data.range(1:data.n_levels) = single(fread(fid,[1, data.n_levels],'single')); % radar altitude layers
    data.T_alt(1:data.T_altcount) = single(fread(fid,[1, data.T_altcount],'single')); % temp prof altitude layers
    data.H_alt(1:data.H_altcount) = single(fread(fid,[1, data.H_altcount],'single')); % hum prof altitude layers
    data.Fr(1:data.n_levels) = int32(fread(fid,[1, data.n_levels],'int')); % range factors
    data.DoppLen = int32(fread(fid,[1, data.no_chirp_seq],'int')); % number of samples in doppler spectra of each chirp sequence
    data.range_offsets = int32(fread(fid,[1, data.no_chirp_seq],'int')) + 1; % chirp sequences start index array in altitude layer array
    data.SeqAvg = int32(fread(fid,[1, data.no_chirp_seq],'int')); % number of averaged chirps within a sequence
    data.SeqIntTime = single(fread(fid,[1, data.no_chirp_seq],'single')); % effective sequence integration time [sec]
    data.dr = single(fread(fid,[1, data.no_chirp_seq],'single')); % chirp sequence range resolution [m]
    data.DoppMax = single(fread(fid,[1, data.no_chirp_seq],'single')); % maximum unambiguious Doppler vel for each chirp sequence
    data.ChanBW = single(fread(fid,[1, data.no_chirp_seq],'single')); % bandwidth of individual radar channel in the sequence [Hz]
    data.ChirpLowIF = int32(fread(fid,[1, data.no_chirp_seq],'int')); % lowest IF frequency in the sequence [Hz]
    data.ChirpHighIF = int32(fread(fid,[1, data.no_chirp_seq],'int')); % highest IF frequency in the sequence [Hz]
    data.RangeMin = int32(fread(fid,[1, data.no_chirp_seq],'int')); % minimum altitude (range) of the sequence [m]
    data.RangeMax = int32(fread(fid,[1, data.no_chirp_seq],'int')); % maximum altitude (range) of the sequence [m]
    data.ChirpFFTSize = int32(fread(fid,[1, data.no_chirp_seq],'int')); % ranging FFT size, must be power of 2
    data.ChirpInvSamples = int32(fread(fid,[1, data.no_chirp_seq],'int')); % number of invalid samples at beginning of chirp
    data.ChirpCenterFr = single(fread(fid,[1, data.no_chirp_seq],'single')); % chirp centre frequency [MHz] at radar transmitter output
    data.ChirpBWFr = single(fread(fid,[1, data.no_chirp_seq],'single')); % chirp bandwidth [MHz] at radar transmitter output
    data.FFTStartInd = int32(fread(fid,[1, data.no_chirp_seq],'int')); % start index of sequence in FFT array
    data.FFTStopInd = int32(fread(fid,[1, data.no_chirp_seq],'int')); % stop index of sequence in FFT array
    data.ChirpFFTNo = int32(fread(fid,[1, data.no_chirp_seq],'int')); % number of FFT range layers in one chirp (usually = 1)
    data.SampRate = int32(fread(fid,1,'int')); % ADC sampling rate [Hz]
    data.MaxRange = int32(fread(fid,1,'int')); % maximum unambiguous range [m]
    data.SupPowLev = int8(fread(fid,1,'char*1')); % flag indicating, if power levelling has been used (0=yes, 1=no)
    data.SpkFilEna = int8(fread(fid,1,'char*1')); % flag indicating, if spike/plankton filter has been used (0=yes, 1=no)
    data.PhaseCorr = int8(fread(fid,1,'char*1')); % flag indicating, if phase correction (dual. pol. radars) has been used (0=yes, 1=no)
    data.RelPowCorr = int8(fread(fid,1,'char*1')); % flag indicating, if relative power correction (dual. pol. radars) has been used (0=yes, 1=no)
    data.FFTWindow = int8(fread(fid,1,'char*1')); % FFT window in use: 0 = SQUARE; 1 = PARZEN; 2 = BLACKMAN; 3 = WELCH; 4 = SLEPIAN2; 5 = SLEPIAN3
    data.FFTInputRng = int32(fread(fid,1,'int')); % ADC input voltage range (+/-)
    data.NoiseFilt = single(fread(fid,1,'single')); % noise filter threshold factor (multiple of STD in Doppler spectra)
    
    % reserved for future use
    fread(fid,25,'int');
    fread(fid,10000,'uint');
    
    data.totsamp = int32(fread(fid,1,'int'));  % total number of samples

    
    % % % % % checks for reasonable header data % % % % % 
    
    % check if binary file contains reasonable number of range bins :
    if any(data.n_levels <= 0 | data.n_levels > 7500 | ~isinteger(data.n_levels)) || isempty(data.n_levels)
        disp(['>>> error opening' infile])
        disp(['>>> file not processed - no of range bins either below 0 or above 7500, empty or not an integer'])
        return
    end

    % check if binary file contains reasonable number of samples :
    if any( data.totsamp <= 0 | data.totsamp > 10000 | ~isinteger(data.totsamp) ) || isempty(data.totsamp)
        disp(['>>> error opening' infile])
        disp(['>>> file not processed - number of samples/profiles below 0 or above 10000, empty or not an integer'])
        return
    end

    % for this file, data transfer did not complete, stopping reading data at last full profile - RG 10.9.2021
    if data.filecode == 889346 && data.starttime == 652622431 && data.endtime == 652625997 ...
        && data.n_levels == 765 && all(data.range_offsets == [1 95 202 388]) ... % strcmp(data.custname, 'Univ. Cologne (Jülich) ')
        && all(data.SeqAvg == [4096 3072 2048 1536]) && all(data.DoppLen == [512 512 512 256])

        data.totsamp = 853;
        
    elseif data.filecode == 889346 && data.starttime == 669751314 && data.endtime == 669754796 ...
        && data.n_levels == 765 && all(data.range_offsets == [1 95 202 388])  ...
        && all(data.SeqAvg == [4096 3072 2048 1536]) && all(data.DoppLen == [512 512 512 256])

        data.totsamp = 86;

    elseif data.filecode == 889346 && data.starttime == 669790822 && data.endtime == 669794389 ...
        && data.n_levels == 765 && all(data.range_offsets == [1 95 202 388])  ...
        && all(data.SeqAvg == [4096 3072 2048 1536]) && all(data.DoppLen == [512 512 512 256])

        data.totsamp = 297;
       
    end
    
    
    % ################################# header ends
       
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% allocate arrays
    
    data.samplen(1:data.totsamp) = int32(0); % sample length in bytes without samplen
    data.time(1:data.totsamp) = int32(0); % time of sample [sec] since 1.1.2001 0:0:0
    data.sampleTms(1:data.totsamp) = int32(0); % milliseconds of sample
    data.QF(1:data.totsamp) = int8(0); % quality flag: bit 4 = ADC saturation, bit 3 = spectral width too high, bit 2 = no transm. power leveling, get bits using dec2bin()
    data.RR(1:data.totsamp) = single(-999); % rain rate [mm/h]
    data.rh(1:data.totsamp) = single(-999); % relative humidity [%]
    data.T_env(1:data.totsamp) = single(-999); % environmental temp [K]
    data.pres(1:data.totsamp) = single(-999); % pressure in [hPa]
    data.ff(1:data.totsamp) = single(-999); % windspeed in [km/h]
    data.fff(1:data.totsamp) = single(-999); % winddirection [°]
    data.u(1:data.totsamp) = single(-999); % voltage [V]
    data.Tb(1:data.totsamp) = single(-999); % brightness temperature [K]
    data.lwp(1:data.totsamp) = single(-999); % liquid water path [g/m³]
    data.powIF(1:data.totsamp) = single(-999); % IF power at ADC [µW]
    data.ele(1:data.totsamp) = single(-999); % eleveation angle [°]
    data.az(1:data.totsamp) = single(-999); % azimuth anlge [°]
    data.status(1:data.totsamp) = single(-999); % status flags: 0/1 heater on/off; 0/10 blower on/off
    data.TransPow(1:data.totsamp) = single(-999); % transmitted power [W]
    data.T_trans(1:data.totsamp) = single(-999); % transmitter temperature [K]
    data.T_rec(1:data.totsamp) = single(-999); % receiver temperature [K]
    data.T_pc(1:data.totsamp) = single(-999); % PC temperature [K]
    data.reserved(1:data.totsamp,1:3) = single(-999);
    data.Tprof(1:data.totsamp,1:data.T_altcount) = single(-999); % temperature profile
    data.Qprof(1:data.totsamp,1:data.H_altcount) = single(-999); % abs hum profile
    data.RHprof(1:data.totsamp,1:data.H_altcount) = single(-999); % rel hum profile
    data.PNv(1:data.totsamp,1:data.n_levels) = single(-999); % total IF power in v-pol measured at the ADC input
    data.SLv(1:data.totsamp,1:data.n_levels) = single(-999); % linear sensitivity limit in Ze units for vertical polarisation
    data.spec(1:data.totsamp,1:data.n_levels,1:max(data.DoppLen)) =  NaN('single'); % vertical pol doppler spectrum linear units % filling with NaNs instead of -999 to avoid having to convert missing values to NaNs later
    if data.DualPol > 0
        data.PNh(1:data.totsamp,1:data.n_levels) = single(-999); % total IF power in h-pol measured at ADT unput
        data.SLh(1:data.totsamp,1:data.n_levels) = single(-999); % linear sensitivity limit in Ze units for horizontal polarisation
        data.spec_hv(1:data.totsamp,1:data.n_levels,1:max(data.DoppLen)) = NaN('single'); % hor pol doppler spectrum linear units
        data.spec_covRe(1:data.totsamp,1:data.n_levels,1:max(data.DoppLen)) = NaN('single'); % real part of covariance spectrum
        data.spec_covIm(1:data.totsamp,1:data.n_levels,1:max(data.DoppLen)) = NaN('single'); % imaginary part of covariance spectrum
    end      
    data.mask(1:data.totsamp,1:data.n_levels) = int8(0); % data.mask array of occupied range cells: 0=not occupied, 1=occupied
    
    if data.CompEna == 2 && data.DualPol > 0
        data.d_spec(1:data.totsamp,1:data.n_levels,1:max(data.DoppLen)) = NaN('single'); % spectral differential reflectivity [dB]
        data.CorrCoeff(1:data.totsamp,1:data.n_levels,1:max(data.DoppLen)) = NaN('single'); % rho_hv, spectral corellation coefficient [0,1]
        data.DiffPh(1:data.totsamp,1:data.n_levels,1:max(data.DoppLen)) = NaN('single'); % spectral differential phase [rad]
    end
    
    if data.DualPol == 2 && data.CompEna == 2
        data.SLDR(1:data.totsamp,1:data.n_levels,1:max(data.DoppLen)) = NaN('single'); % compressed spectral slanted LDR [dB]
        data.SCorrCoeff(1:data.totsamp,1:data.n_levels,1:max(data.DoppLen)) = NaN('single'); % compressed spectral slanted correlation coefficient [0,1]
        if data.CompEna == 2
            data.KDP(1:data.totsamp,1:data.n_levels) = single(-999); % specific differential phase shift [rad/km]
            data.DiffAtt(1:data.totsamp,1:data.n_levels) = single(-999); % differential attenuation [dB/km]             
        end % data.CompEna        
    end % data.DualPol
    
    if data.CompEna > 0
        data.VNoisePow_mean(1:data.totsamp,1:data.n_levels) = single(-999); % integrated Doppler spectrum noise power in v-pol [Ze]
        if data.DualPol > 0
            data.HNoisePow_mean(1:data.totsamp,1:data.n_levels) = single(-999); % integrated Doppler spectrum noise power in h-pol [Ze]
        end
    end
                    
    if data.AntiAlias == 1
        data.Aliasmask(1:data.totsamp,1:data.n_levels) = int8(0); % data.mask array if aliasing has been applied: 0=not apllied, 1=apllied
        data.MinVel(1:data.totsamp,1:data.n_levels) = single(-999.); % minimum velocity in Doppler spectrum [m/s]
    end
    
    
    %############################## end allocating
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start reading the rest of the file
    for i = 1:data.totsamp
        
        data.samplen(i) = fread(fid,1,'int');
        data.time(i) = fread(fid,1,'uint');
        data.sampleTms(i) = fread(fid,1,'int');
        data.QF(i) = int8(fread(fid,1,'char*1'));
        data.RR(i) = fread(fid,1,'single');
        data.rh(i) = fread(fid,1,'single');
        data.T_env(i) = fread(fid,1,'single');
        data.pres(i) = fread(fid,1,'single');
        data.ff(i) = fread(fid,1,'single');
        data.fff(i) = fread(fid,1,'single');
        data.u(i) = fread(fid,1,'single');
        data.Tb(i) = fread(fid,1,'single');
        data.lwp(i) = fread(fid,1,'single');
        data.powIF(i) = fread(fid,1,'single');
        data.ele(i) = fread(fid,1,'single');
        data.az(i) = fread(fid,1,'single');
        data.status(i) = fread(fid,1,'single');
        data.TransPow(i) = fread(fid,1,'single');
        data.T_trans(i) = fread(fid,1,'single');
        data.T_rec(i) = fread(fid,1,'single');
        data.T_pc(i) = fread(fid,1,'single');
        % reserved: 
        data.reserved(i,1:3) = fread(fid,3,'single'); % reserved
        %fread(fid,3,'single');
        
        data.Tprof(i,1:data.T_altcount) = fread(fid,data.T_altcount,'single'); % temp prof
        data.Qprof(i,1:data.H_altcount) = fread(fid,data.H_altcount,'single'); % abs hum prof
        data.RHprof(i,1:data.H_altcount) = fread(fid,data.H_altcount,'single'); % rel hum prof

        data.PNv(i,1:data.n_levels) = fread(fid,[1, data.n_levels],'single');
        if data.DualPol > 0
            data.PNh(i,1:data.n_levels) = fread(fid,[1, data.n_levels],'single');
        end
            
            
        try
            data.SLv(i,1:data.n_levels) = fread(fid,[1, data.n_levels],'single');
            if data.DualPol > 0
                data.SLh(i,1:data.n_levels) = fread(fid,[1, data.n_levels],'single');
            end

        catch
             if feof(fid)
                endoffileerrormessage
                return
             end
        end

            
        
        data.mask(i,1:data.n_levels) = int8(fread(fid,[1, data.n_levels],'char*1'));        
        
        for j = 1:data.n_levels
                        
            if data.mask(i,j) == 1
                

                
                fread(fid,1,'int'); % number of bytes of the followng spectral block
                
                % spectra
                if data.CompEna == 0 

                    % find number of current chirp index
                    chirp_idx = int32(find(data.range_offsets(2:end) - j > 0,1,'first'));
                    if isempty(chirp_idx) % then j is within last chirp
                        chirp_idx = data.no_chirp_seq;
                    end

                     try

                        data.spec(i,j,1:data.DoppLen(chirp_idx)) = fread(fid,[1,data.DoppLen(chirp_idx)],'single'); % spectra
                        if data.DualPol > 0
                            data.spec_hv(i,j,1:data.DoppLen(chirp_idx)) = fread(fid,[1,data.DoppLen(chirp_idx)],'single'); % spectra 
                            data.spec_covRe(i,j,1:data.DoppLen(chirp_idx)) = fread(fid,[1,data.DoppLen(chirp_idx)],'single'); % spectra 
                            data.spec_covIm(i,j,1:data.DoppLen(chirp_idx)) = fread(fid,[1,data.DoppLen(chirp_idx)],'single'); % spectra 
                        end

                    catch
                         if feof(fid)
                            endoffileerrormessage
                            return
                         end
                    end
                    
                
                else %  data.CompEna > 0 
                    
                    Nblocks = int8(fread(fid,1,'char*1')); % number of blocks in spectra
                    MinBkIdx = fread(fid,[1, Nblocks],'int16') + 1; % minimum indexes of blocks
                    MaxBkIdx = fread(fid,[1, Nblocks],'int16') + 1; % maximum indexes of blocks
                    
                    for jj = 1:Nblocks
                        try
                            data.spec(i,j,MinBkIdx(jj):MaxBkIdx(jj)) = fread(fid,[1,MaxBkIdx(jj)-MinBkIdx(jj)+1],'single');
                        catch
                             if feof(fid)
                                endoffileerrormessage
                                return
                             end
                        end
                    end % jj
                    
                    if data.DualPol > 0
                         
                        try

                            for jj = 1:Nblocks
                                data.spec_hv(i,j,MinBkIdx(jj):MaxBkIdx(jj)) = fread(fid,[1,MaxBkIdx(jj)-MinBkIdx(jj)+1],'single');
                            end
                        
                            for jj = 1:Nblocks
                                data.spec_covRe(i,j,MinBkIdx(jj):MaxBkIdx(jj)) = fread(fid,[1,MaxBkIdx(jj)-MinBkIdx(jj)+1],'single');
                            end
                        
                            for jj = 1:Nblocks
                                data.spec_covIm(i,j,MinBkIdx(jj):MaxBkIdx(jj)) = fread(fid,[1,MaxBkIdx(jj)-MinBkIdx(jj)+1],'single');
                            end
                            
                        catch
                            disp('what?')
                            if feof(fid)
                                endoffileerrormessage
                                return
                            end
                        end
                    end % if
                        
                    if data.CompEna == 2
                        
                        disp('!!WARNING!! in read_lv0_v3: the combination of DualPol > 0 and CompEna == 2 has not been tested')

                        for jj = 1:Nblocks
                            data.d_spec(i,j,MinBkIdx(jj):MaxBkIdx(jj)) = fread(fid,[1,MaxBkIdx(jj)-MinBkIdx(jj)+1],'single');
                        end
                        for jj = 1:Nblocks
                            data.CorrCoeff(i,j,MinBkIdx(jj):MaxBkIdx(jj)) = fread(fid,[1,MaxBkIdx(jj)-MinBkIdx(jj)+1],'single');
                        end
                        for jj = 1:Nblocks
                            data.DiffPh(i,j,MinBkIdx(jj):MaxBkIdx(jj)) = fread(fid,[1,MaxBkIdx(jj)-MinBkIdx(jj)+1],'single');
                        end

                        if data.DualPol == 2
                            for jj = 1:Nblocks
                                data.SLDR(i,j,MinBkIdx(jj):MaxBkIdx(jj)) = fread(fid,[1,MaxBkIdx(jj)-MinBkIdx(jj)+1],'single');
                            end
                            for jj = 1:Nblocks
                                data.SCorrCoeff(i,j,MinBkIdx(jj):MaxBkIdx(jj)) = fread(fid,[1,MaxBkIdx(jj)-MinBkIdx(jj)+1],'single');
                            end 
                        end
                           
                    end %fi
                    
                    if data.DualPol == 2 && data.CompEna == 2
                        data.KDP(i,j) = fread(fid,1,'single');
                        data.DiffAtt(i,j) = fread(fid,1,'single');
                    end
                    
                    data.VNoisePow_mean(i,j) = fread(fid,1,'single');
                    if data.DualPol > 0
                        data.HNoisePow_mean(i,j) = fread(fid,1,'single');
                    end
                    
                    if data.AntiAlias == 1
                        data.Aliasmask(i,j) = int8(fread(fid,1,'char*1'));
                        data.MinVel(i,j) = fread(fid,1,'single');
                    end
                                        
                end % data.CompEna == 0

            
            end % if data.mask(i,j)
        end % j
        
    end % i
    
    fclose(fid);                    
          
    data.readerror = false; % succesfull in reading file!                 
    

end %function


function endoffileerrormessage()

disp(' !!! WARNING !!!')
disp('    -> End of file reached unexpectedly. It is possible that lv0 file is corrupted. Check!')
end
