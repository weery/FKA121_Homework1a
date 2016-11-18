% Plot the average temperature % close previous figure

clear all

data = importdata('data.dat');
data2= importdata('data2.dat');


% plot
figure;
startx = 20;
plot(data(startx:end,1),data(startx:end,2).*data(startx:end,1).^2,'-')

title('Static structure factor')
xlabel('q')
ylabel('S(q)')