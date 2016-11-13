
% plot the displacements
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('velocities1.dat');

sizeD=size(data);

%plot
figure;
hold on
for t=1:1:1000
    for i=1:256
        hold on
        plot3(data(t,1+(i-1)*3),data(t,2+(i-1)*3),data(t,3+(i-1)*3),'X')
    end
    pause(0.1)
    clf
end
% labels
xlabel('Time / [dim. unit]');
ylabel('Displacement / [dim. unit]');

% legend
legend('Atom 1','Atom 2','Atom 3');

% axis limits

