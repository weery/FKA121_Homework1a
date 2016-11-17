% Plot the average temperature
close % close previous figure

data = importdata('data.dat');


% plot
figure;
%%plot(data(:,2))



plot(data(:,1),data(:,2))
