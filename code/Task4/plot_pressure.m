% Plot the average temperature
close % close previous figure
clear all

data = importdata('pressure.dat');

% plot
figure;
plot(data(1:12000,1), data(1:12000,2), '-')

hold on
target_temperature = 101325e-11/1.602;
plot(data(1:12000,1),target_temperature*ones(1,12000),'--')

% legend
legend({'Instantanous pressure','Target pressure'},'interpreter','latex','location','northeast')
% 
% labels
xlabel('Time / [ps]','interpreter','latex')
ylabel('Pressure / [eV/\AA]','interpreter','latex')
% title
title('Pressure equilibration','interpreter','latex')