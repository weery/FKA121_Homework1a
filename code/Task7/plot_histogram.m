% Plot the average temperature
close % close previous figure
clear all
data = importdata('histogram.dat');


% plot
figure;
%%plot(data(:,2))



g=data(:,2)./data(:,3);


plot(data(:,1),g)

xlim([0,data(end,1)])
title('Radial distribution function g(r)')
xlabel('r [Å]')
ylabel('g(r)')