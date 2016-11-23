% Plot the average temperature
close % close previous figure
clear all
data = importdata('histogram.dat');


% plot
figure;
%%plot(data(:,2))


% Rätt faktor? antal tidssteg och antal atomer
g=data(:,3)./data(:,4)/1000/256;


plot(data(:,1),g)

x=57;
summation =sum(g(1:x)*4*pi.*data(1:x,1).^2)/256;


xlim([0,data(end,1)])
title('Radial distribution function g(r)')
xlabel('$r$ / [\r{A}]', 'Interpreter', 'latex')
ylabel('$g(r)$', 'Interpreter', 'latex')