close % plot the energies
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('energy.dat');
k_b=0.00008617;
N=256;
%plot 
figure;
const=mean((data(:,3)-mean(data(:,3))).^2);
c_v=(3*N*k_b/2)/(1-2/(3*N*k_b^2*500^2)*const);


% labels
xlabel('Time / [ASU]');
ylabel('Energy / [ASU]');

% legend
legend('Total energy','Potential energy','Kinetic energy');

title('Awesome title')
