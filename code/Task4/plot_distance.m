% plot the displacements
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('displacement.dat');

sizeD=size(data);

%plot
figure;
hold on
dist = zeros(10000,1);
for j=1:256
    tra=[data(1,1+(j-1)*3)-data(:,1+(j-1)*3),data(:,2+(j-1)*3)-data(:,2+(j-1)*3),data(:,3+(j-1)*3)-data(:,3+(j-1)*3)];

    for i=1:10000
        dist(i)=dist(i)+norm(tra(i,:))/4.046;
    end
end
    plot(1:10000,dist,'-')