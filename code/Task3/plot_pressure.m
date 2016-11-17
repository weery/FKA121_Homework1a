% Plot pressure
close % close previous figure

data = importdata('pressure_avg.dat');

% plot
figure;
plot(data(4:end,1), data(4:end,2), '-');

% labels
xlabel('Time / [ps]')
ylabel('Pressure / [eV/A]')
% title
title('Pressure average')