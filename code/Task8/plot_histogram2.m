% Plot the average temperature
close % close previous figure
clear all
data = importdata('../Task7/histogram.dat');


% plot
figure;
%%plot(data(:,2))

g=data(:,3)./data(:,4)/1000/32;

%plot(data(:,1),g)

nDats=25;

integr=@(q)sum(data(1:nDats,1).^2.*(g(1:nDats)-1).*sin(q*data(1:nDats,1))./(data(1:nDats,1)*q))/((data(nDats,1)-data(1,1))*nDats);

s=@(q)1+4*pi*100*integr(q);


c=linspace(5,25,5000);

for i=1:5000
    b(i)=s(c(i));
end
plot(c,b)

xlim([c(1), c(end)])

title('Static structure factor')
xlabel('q')
ylabel('S(q)')