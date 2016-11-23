% Plot the average temperature
close % close previous figure
clear all
data = importdata('histogram.dat');


% plot
figure;
%%plot(data(:,2))


% Rätt faktor? antal tidssteg och antal atomer
g=data(:,3)./data(:,4)/1000/32;


plot(data(:,1),g)

xlim([0,data(end,1)])
title('Radial distribution function g(r)')
xlabel('$r$ / [\r{A}]', 'Interpreter', 'latex')
ylabel('$g(r)$', 'Interpreter', 'latex')