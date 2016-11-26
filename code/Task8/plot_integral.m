% Plot the average temperature
 % close previous figure
clear all
data = importdata('../Task7/histogram.dat');

n=256;
% Divide by number of atoms and timesteps
g=data(:,3)./data(:,4)/10000/n;

%plot(data(:,1),g)

% about half the cell length
nDats=125;

rad=data(1:nDats,1);

integr=@(q)sum(rad.^2.*(g(1:nDats)-1).*sin(q*rad)./(rad*q));

s=@(q)1+4*pi*1/n*integr(q);

q=linspace(0,25,5000);

for i=1:5000
    val(i)=s(q(i));
end
plot(q,val)

xlim([q(1), q(end)])

title('Integrated static structure factor','interpreter','latex')
xlabel('$q$','interpreter','latex')
ylabel('$s(q)$','interpreter','latex')