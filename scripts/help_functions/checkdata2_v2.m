function data = checkdata2_v2(data)
% Sometimes issues with time stamp due to GPS failures (presumably). 
% Checking if time goes "backwards" or duplicate time stamps presents, and 
% correct time stamps where necessary. Lukas had the great idea of checking 
% if we can find a corresponding forward jump and correct the time stamp.
% RG 1.12.2022     

data.timeshift = zeros(size(data.time)); % record how much time shifted 
% positive value means that the amount was added to the original time
% stamp, to get back the original timestamp, calculate 
% timestamp(=time+sampleTms*1e-3) - timeshift

if ~any(diff(double(data.time) + double(data.sampleTms).*1e-3) <= 0 ) 
    % no problem found
    return
end

timestamp = double(data.time) + double(data.sampleTms).*1e-3;

% find temporal resolution of measurement: assume most common difference
% between time stamps is the expected duration of measuring each profile
time_res = median(diff(timestamp));

%% correct backward jumps 

% one file (joyrad94_20170522150002_P05_ZEN.lv0) found that has problems 
% to the extend that the algorithm cannot handle them: creating an
% exception to deal with this 

if data.time(1) == 517158002 && data.n_levels == 1021 && all(data.range_offsets == [1 76 226 494])

    data = special_case(data, timestamp, 565:625);

% same issue for joyrad94_20170527210000_P05_ZEN.lv0
elseif data.time(1) == 517611600 && data.n_levels == 1021 && all(data.range_offsets == [1 76 226 494])

    data = special_case(data, timestamp, 225:276);

else 

    % also special treatment needed in joyrad94_20170623180001_P05_ZEN.lv0    
	if data.time(1) == 519933601 && data.n_levels == 1021 && all(data.range_offsets == [1 76 226 494])    
        data = special_case(data, timestamp, 1224:1234);
        timestamp = double(data.time) + double(data.sampleTms).*1e-3;
        
    elseif data.time(1) == 519901200 && data.n_levels == 1021 && all(data.range_offsets == [1 76 226 494])    
        data = special_case(data, timestamp, 19:23);
        timestamp = double(data.time) + double(data.sampleTms).*1e-3;
        
    % same issue for joyrad94_20170719130002_P05_ZEN.lv0
    elseif data.time(1) == 522162002 && data.n_levels == 1021 && all(data.range_offsets == [1 76 226 494])
        data = special_case(data, timestamp, 354:414);
        timestamp = double(data.time) + double(data.sampleTms).*1e-3;
        
    end
    
    % within the function correct_backward_jump, if there are jumps very close
    % to each other only one of them is evaluated, since experience showed that
    % they are often connected and it is enough to correct the last one. 
    % For the rare cases where the problems is not fixed like this, the
    % function is called in a loop, so that the previously ignored jumps also
    % get corrected.

    nn = 0;
    while any(diff(timestamp) < 0 )
        nn = nn + 1;

        data = correct_backward_jump(data, timestamp, time_res);

        % time most likely edited in function, update time stamp
        timestamp = double(data.time) + double(data.sampleTms).*1e-3;    

        if nn > 5
            disp('backward jump in time not resolved')
            break
        end

    end

end


% check that no data point outside of the hour being processed
data = checkdata1(data);

% if problem solved, return
if ~any(diff(double(data.time) + double(data.sampleTms).*1e-3) <= 0 )
    return
end


%% for dublicate time stamp: move the second of the duplicates to be 
% half way between previous and next one

timestamp = double(data.time) + double(data.sampleTms).*1e-3;
ind = find(diff(timestamp) == 0 )+1; % add +1 to index the second occurrence of the duplicate

for tt = ind

    if tt == length(data.time) % duplicate time stamp at the end of the file
        add_time = time_res;
        
    else
        % calculate mean of duplicate and next time stamp after dubplicate stamp
        add_time = (timestamp(tt+1) - timestamp(tt-1)) / 2;
        add_time = min(add_time, time_res); % don't add more than the temporal resolution
    end
    
    new_time = timestamp(tt) + add_time;
    
    data.time(tt) = floor(new_time);
    data.sampleTms(tt) = (new_time - floor(new_time)).*1e3;

    data.timeshift(tt) = data.timeshift(tt) + add_time;
        
    data.totsampchangelabel = 1; % note for later that issues with time found

end


% check that no data point outside of the hour being processed
% probably not necessary here, but better safe than sorry
data = checkdata1(data);


if ~any(diff(double(data.time) + double(data.sampleTms).*1e-3) <= 0 ) 
    return
end


%% still a problem with time? 
    
disp('WARNING ! This not programmed yet')

end % function


function data = correct_backward_jump(data, timestamp, time_res)


    % find backward jumps
    ind = find(diff(timestamp) < 0 );

    % in a couple of files two jumps occur close to each other, which messes up
    % with the algorithm. solution in the cases found was to ignore first jump
    if any(diff(ind) < 5)
        ind( diff(ind) < 5 ) = [];
    end

    for tt = ind

        figure, plot(timestamp, '.'), xlim([tt-20, tt+10]), hold on

        % 1. check if backward jump associated with a corresponding forward jump before

        % starting from the time stamp after the backward jump, calculate
        % backwards what the time stamps would be base corresponding to the
        % temporal resolution, and compare to the logged timestamp

        testtime = timestamp(tt+1) - [1:tt].*time_res;

        %  where time stamp within +- 0.3 sec from testtime
        ix_close = abs(testtime - flip(timestamp(1:tt))) < 0.3;    
        ix = find(ix_close, 1); % find first ocurrence

        % found possible match going backwards:
        % criteria 1: a match found in the past 
        % criteria 2: at the point when timestamps start to align, is there a forward jump of at least time_res+0.3
        %   second criteria is to avoid drifting towards similar time due to
        %   inprecision of the time resolution
        %
        % -> consider that found the corresponding forward jump -

        if ~isempty(ix) && (timestamp(tt-ix+2)-timestamp(tt-ix+1)) > time_res+0.3

            new_time =  flip(testtime(1:ix-1)); % time to be filled in

            ix_ed = (tt-ix+2:tt); % index of original time array to be edited

            data.timeshift(ix_ed) =  data.timeshift(ix_ed) + ... % calculate cumulatively, since possible that same time stamp gets modified more than once
                new_time - ( double(data.time(ix_ed)) + double(data.sampleTms(ix_ed)).*1e-3 );

            data.time(ix_ed) = floor(new_time);
            data.sampleTms(ix_ed) = (new_time-floor(new_time) ) .*1e3;

            data.totsampchangelabel = 1; % note for later that issues with time found

        % 2. if no corresponding forward jump before, correct forward
        else 

            %  where time stamp within +- 0.3 sec from testtime
            testtime = timestamp(tt) + ([tt : length(data.time)-1] -tt+1).*time_res;

            ix_close = abs(testtime - timestamp(tt+1:end)) < 0.3;
            ix = find(ix_close, 1); % find first ocurrence

            % found possible match going forwards
            % as above, also check that a associated with a jump
            if any(ix_close) && (timestamp(tt+ix)-timestamp(tt+ix-1)) > time_res+0.3

                new_time = testtime( 1: ix-1); % time to be filled in

                ix_ed = ( tt+1: tt+ix-1); % index of original time array to be edited

                data.timeshift(ix_ed) =  data.timeshift(ix_ed) + ... % calculate cumulatively, since possible that same time stamp gets modified more than once
                    new_time - ( double(data.time(ix_ed)) + double(data.sampleTms(ix_ed)).*1e-3 );

                data.time(ix_ed) = floor(new_time);
                data.sampleTms(ix_ed) = (new_time-floor(new_time) ) .*1e3;

                data.totsampchangelabel = 1; % note for later that issues with time found



            else % no matching time found afterwards (and also not before, because that checked above)
                tix = tt; 

                % now it is of course possible that now the backwards jump has
                % moved

                % continue until end of file or until a forward jump occurs 
                % (could be because of the internal calibration)
                while (timestamp(tix+1) - timestamp(tix)) < 0 

                    % shift time stamp to be the time stamp before jump + time res
                    new_time = timestamp(tix) + time_res;

                    data.timeshift(tix+1) = data.timeshift(tix+1) + new_time - ...
                         ( double(data.time(tix+1)) + double(data.sampleTms(tix+1)).*1e-3 );


                    % need to modify also variable timestamp, as this used in
                    % the statement of the while loop
                    data.time(tix+1) = floor(new_time);
                    data.sampleTms(tix+1) = (new_time-floor(new_time) ) .*1e3;
                    timestamp(tix+1) =  new_time;

                    data.totsampchangelabel = 1; % note for later that issues with time found

                    tix = tix+1; % move to check next time stamp
                    
                    if tix == length(data.time) % reached end of file
                        break
                    end
                end

            end

        end % fi

        flag_edittime = data.timeshift ~= 0 ;

        plot(find(flag_edittime),  double(data.time(flag_edittime)) + double(data.sampleTms(flag_edittime)).*1e-3 , 'or')
        legend('original timestamp', 'modified time stamp', 'location', 'northwest')
        xlim([tt-20, find(flag_edittime, 1, 'last')+10])

        this_date = datestr(double(data.time(1))/3600/24 + datenum([2001,1,1,0,0,0]), 'yyyymmddHH');

        title(this_date)
        savefig(['timetstamp_shiftcheck_' this_date '_' num2str(tt)])
        close

    end % tt


end % function


function data = special_case(data, timestamp, ind_out)

    ind = false(size(timestamp));
    ind(ind_out) = 1;

    figure, plot(timestamp, '.'), xlim([find(ind,1)-20, find(ind,1,'last')+20]), hold on

    new_time = interp1(find(~ind), timestamp(~ind), find(ind), 'linear');
    
    data.timeshift(ind) =  data.timeshift(ind) + ... % calculate cumulatively, since possible that same time stamp gets modified more than once
        new_time - ( double(data.time(ind)) + double(data.sampleTms(ind)).*1e-3 );

    data.time(ind) = floor(new_time);
    data.sampleTms(ind) = (new_time-floor(new_time) ) .*1e3;

    data.totsampchangelabel = 1; % note for later that issues with time found

    
    flag_edittime = data.timeshift ~= 0 ;

    plot(find(flag_edittime),  double(data.time(flag_edittime)) + double(data.sampleTms(flag_edittime)).*1e-3 , 'or')
    legend('original timestamp', 'modified time stamp', 'location', 'northwest')
        
    this_date = datestr(double(data.time(1))/3600/24 + datenum([2001,1,1,0,0,0]), 'yyyymmddHH');

    title(this_date)
    savefig(['timetstamp_shiftcheck_' this_date '_specialcase'])
    close

end
