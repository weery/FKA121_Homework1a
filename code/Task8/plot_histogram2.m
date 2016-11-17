% Plot the average temperature
close % close previous figure
clear all
data = importdata('histogram.dat');


% plot
figure;
%%plot(data(:,2))

g=data(:,3)./data(:,4);

%plot(data(:,1),g)

nDats=20;

integr=@(q)sum(data(1:nDats,1).^2.*(g(1:nDats)-1).*sin(q*data(1:nDats,1))./(data(1:nDats,1)*q))/((data(nDats,1)-data(1,1))*nDats);

s=@(q)1+4*pi*100*integr(q);



for i=1:500
    b(i)=s(i/150);
end
c=linspace(0,5,500);
plot(c,b)

title('Static structure factor')
xlabel('q')
ylabel('S(q)')