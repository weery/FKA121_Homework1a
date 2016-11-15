close % plot the energies
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('energy.dat');

%plot 
figure;
plot(data(:,1:20:end),data(:,2:20:end)+data(:,3:20:end),'-')
hold on
plot(data(:,1:20:end),data(:,2:20:end),'-')
plot(data(:,1:20:end),data(:,3:20:end),'-')

% labels
xlabel('Time / [ASU]');
ylabel('Energy / [ASU]');

% legend
legend('Total energy','Potential energy','Kinetic energy');

title('Awesome title')
