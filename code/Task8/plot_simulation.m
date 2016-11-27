% Plot the average temperature % close previous figure

clear all

data = importdata('sq.dat');
data2= importdata('sq_bin.dat');

% plot
%figure;
x_start = 265722;
width=50;
stuff = (data(x_start:(x_start+width),1)+data(x_start:(x_start+width),2));
stuff2 = (data(x_start:(x_start+width),3));

%stuff(x_start)=0;
%histogram(stuff,1000)


zeroeses=find (data(:,3)==0);
data(zeroeses,:)=[];

[a,b]=sort(data(:,3));

plot(a,data(b,1)+data(b,2),'.')


title('Simulated Static structure factor','interpreter','latex')
xlabel('q','interpreter','latex')
ylabel('S(q)','interpreter','latex')