% Plot pressure
close % close previous figure

data = importdata('pressure.dat');

% plot
figure;
plot(data(4:end,1), data(4:end,2), '-');

% labels
xlabel('Time / [ps]')
ylabel('Average pressure / [eV/A]')
% title
title('Average pressure')