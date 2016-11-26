% Plot the average temperature
close % close previous figure
clear all

data = importdata('temperature.dat');

% plot
figure;
plot(data(1:8000,1), data(1:8000,2), '-')

hold on
target_temperature = 500+273.15;
plot(data(1:8000,1),target_temperature*ones(1,8000),'--')

% legend
legend({'Instantanous temperature','Target temperature'},'interpreter','latex','location','east')
% 
% labels
xlabel('Time / [ps]','interpreter','latex')
ylabel('Temperature / [K]','interpreter','latex')
% title
title('Temperature equilibration','interpreter','latex')