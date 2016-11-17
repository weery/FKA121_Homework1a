% Plot the average temperature
close % close previous figure

data = importdata('temperature_avg.dat');

% plot
figure;
plot(data(:,1), data(:,2), '-')

% labels
xlabel('Time / [ps]')
ylabel('Temperature / [K]')
% title
title('Temperature average')