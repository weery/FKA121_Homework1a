%clear all, close all, clc
clf
addpath( genpath('.') );

s = strcat('energy0.00',int2str(5),'.dat');
    data = importdata(s);
    hold on
    plot(data(:,1),data(:,3));

for i = 2:8
    s = strcat('energy0.0',int2str(5*i),'.dat');
    data = importdata(s);
    hold on
    plot(data(:,1),data(:,3));
end


text(5.2,28000000,'0.03 \uparrow')

text(1,5000000,'\downarrow 0.005')
text(2.5,5000000,'\downarrow 0.01')
text(3.8,5000000,'\downarrow 0.015')
text(5,5000000,'\downarrow 0.02')
text(6.5,5000000,'\downarrow 0.025')

text(3,32000000,'0.035 \downarrow')
text(4,40000000,'0.04 \downarrow')
% labels
xlabel('Time / [ps]','Interpreter','Latex');
ylabel('Energy / [eV]','Interpreter','Latex');

title('Total energy over time','Interpreter','Latex')
