clear all, close all, clc

data=importdata('sq.dat');


len=size(data);
data(1,2)=0;
[Y,I]=sort(data(:,2));

dataPlot=data(I,1);



plot(dataPlot)



