% Plot the average temperature % close previous figure

clear all

data = importdata('data.dat');
data2= importdata('data2.dat');


% plot
figure;
plot(data(2:10,1),data(2:10,2),'-')

title('Static structure factor')
xlabel('q')
ylabel('S(q)')