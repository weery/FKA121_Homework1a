% Plot the average temperature % close previous figure

clear all

data = importdata('data.dat');
data2= importdata('data2.dat');


% plot
figure;
startx = 1;
plot(data(startx:end,1),data(startx:end,2),'-')

title('Simulated Static structure factor')
xlabel('q')
ylabel('S(q)')