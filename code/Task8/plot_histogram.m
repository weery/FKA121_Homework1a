% Plot the average temperature
close % close previous figure

data = importdata('histogram.dat');


% plot
figure;
%%plot(data(:,2))

plot(data(:,1),data(:,3)./data(:,4)/256)

xlim([0,data(end,1)])