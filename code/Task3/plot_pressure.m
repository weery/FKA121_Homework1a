edi% Plot pressure
close % close previous figure

data = importdata('pressure.dat');

% plot
figure;
plot(data(:,1), data(:,2), '-');

% labels
xlabel('Time / [ps]')
ylabel('Average pressure / [eV/A]')
% title
title('Average pressure')