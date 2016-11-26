% Plot the average temperature
 % close previous figure
clear all
data = importdata('../Task7/histogram.dat');

n=256;
% Divide by number of atoms and timesteps
g=data(:,3)./data(:,4)/10000/n;

%plot(data(:,1),g)

nDats=125;


rad=data(1:nDats,1);

integr=@(q)sum(rad.^2.*(g(1:nDats)-1).*sin(q*rad)./(rad*q));
%)/((data(nDats,1)-data(1,1))*nDats);

s=@(q)1+4*pi*1/n*integr(q);


c=linspace(0,25,5000);

for i=1:5000
    b(i)=s(c(i));
end
plot(c,b)

xlim([c(1), c(end)])

title('Static structure factor')
xlabel('q')
ylabel('S(q)')