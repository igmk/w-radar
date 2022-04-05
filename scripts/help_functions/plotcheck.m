function plotcheck(data)

radar_time = datenum([2001 1 1 0 0 0]) + double(data.time)/(86400);

figure('Units', 'pixels', 'Position', [0 0 800 1000]);

% Ze
subplot(3,1,1)
TimeHeightPlot(radar_time, data.range, 10.*log10( data.Ze'), 'Z_e [dBz]')

% Vm
subplot(3,1,2)
TimeHeightPlot(radar_time, data.range, data.vm', 'V_m [m/s]')

% spectrum width
subplot(3,1,3)
TimeHeightPlot(radar_time, data.range, data.sigma', 'Sepctrum width [m/s]')

savefig(['CheckTimeHeight' datestr(radar_time(1), 'yyyymmddHH') ] )
close

function TimeHeightPlot(tim, r, dat, collbl)

pcolor(tim, r, dat)
shading flat
datetick

cb = colorbar;
ylabel(cb, collbl)
shading flat
ylabel('Range [m]');
xlabel('Time');