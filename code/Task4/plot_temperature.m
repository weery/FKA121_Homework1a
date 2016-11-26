% Plot the average temperature
close % close previous figure
clear all

data = importdata('temperature.dat');

% plot
figure;
plot(data(1:12000,1), data(1:12000,2), '-')

hold on
target_temperature1 = 1500+273.15;
target_temperature2 = 700+273.15;
plot(data(1:6000,1),target_temperature1*ones(1,6000),'-')
plot(data(6001:12000,1),target_temperature2*ones(1,6000),'-')

% legend
legend({'Instantanous temperature','Target temperature melting phase', 'Target temperature'},'interpreter','latex','location','southeast')
% 
% labels
xlabel('Time / [ps]','interpreter','latex')
ylabel('Temperature / [K]','interpreter','latex')
% title
title('Temperature equilibration','interpreter','latex')