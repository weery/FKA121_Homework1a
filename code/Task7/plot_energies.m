close % plot the energies
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('energy.dat');

%plot 
figure;
plot(data(:,1),data(:,2)+data(:,3),'-')
hold on
plot(data(:,1),data(:,2),'-')
plot(data(:,1),data(:,3),'-')

% labels
xlabel('Time / [ASU]');
ylabel('Energy / [ASU]');

% legend
legend('Total energy','Potential energy','Kinetic energy');

title('Awesome title')