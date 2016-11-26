% Plot the average temperature
close % close previous figure
clear all
data = importdata('histogram.dat');
% plot
figure;
% Rätt faktor antal tidssteg och antal atomer
g=data(:,3)./data(:,4)/10000/256;

plot(data(:,1),g)

%Integrate to first minimum
first_minimum_idx=55;
summation =sum(g(2:first_minimum_idx)*4*pi.*data(2:first_minimum_idx,1).^2)/256;

xlim([0,data(end,1)])
title('Radial distribution function $g(r)$','interpreter','latex')
xlabel('$r$ / [\r{A}]', 'Interpreter', 'latex')
ylabel('$g(r)$', 'Interpreter', 'latex')