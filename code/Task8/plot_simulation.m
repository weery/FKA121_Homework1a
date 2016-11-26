% Plot the average temperature % close previous figure

clear all

data = importdata('sq.dat');
data2= importdata('sq_bin.dat');


% plot
figure;
startx = 1;
plot(data(startx:end,1),data(startx:end,2),'-')


title('Simulated Static structure factor','interpreter','latex')
xlabel('q','interpreter','latex')
ylabel('S(q)','interpreter','latex')