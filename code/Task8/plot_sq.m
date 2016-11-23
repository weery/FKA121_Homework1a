clear all, close all, clc

data=importdata('sq.dat');

[Y,I]=sort(data(:,2));

dataPlot=data(I,1);

binSize=0.4;
nBins= ceil(max(Y)/binSize);

bins=zeros(1,nBins);

for j=1:length(Y)
    len=Y(j);
    for i=1:nBins
        if len< i*binSize
            bins(i)=+dataPlot(i);
            break;
        end
    end    
end

plot(bins)