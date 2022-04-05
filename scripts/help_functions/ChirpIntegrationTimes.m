function inttimes = ChirpIntegrationTimes(Vnyq, f0, N_c)
% chirp integration time, i.e. the time required for each chirp, is
% calculated from chirp parameters
% Rosa Gierens 4.12.2020

c = 3e8; % speed of light

% chirp repetition frequency
F_c = 4 * Vnyq * f0*1e9 / c;


% chirp duration T_c = 1 / F_c

% chirp sequency duration T_int = T_n * N_c
inttimes =  1 ./ F_c .* double(N_c);


%data.ChirpIntTime 