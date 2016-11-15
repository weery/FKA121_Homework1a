% Plot the average temperature
close % close previous figure

data = importdata('temperature.dat');

% plot
figure;
plot(data(:,1), data(:,2), '-')

% labels
xlabel('Time / [ps]');
ylabel('Average temperature / [K]');