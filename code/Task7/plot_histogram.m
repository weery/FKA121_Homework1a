% Plot the average temperature
close % close previous figure
clear all
data = importdata('histogram.dat');


% plot
figure;
%%plot(data(:,2))


% R�tt faktor? antal tidssteg och antal atomer
g=data(:,3)./data(:,4)/1000/256;


plot(data(:,1),g)

[rad,top]=max(g);

summation =(sum(g(1:top))+sum(g(1:top+1)))/2;

xlim([0,data(end,1)])
title('Radial distribution function g(r)')
xlabel('$r$ / [\r{A}]', 'Interpreter', 'latex')
ylabel('$g(r)$', 'Interpreter', 'latex')