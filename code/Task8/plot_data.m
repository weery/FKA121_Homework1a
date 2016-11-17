% Plot the average temperature
close % close previous figure
close all
clear all

data = importdata('data.dat');
data2= importdata('data2.dat');


% plot
figure;




plot(data2(:,2),data2(:,1),'-')