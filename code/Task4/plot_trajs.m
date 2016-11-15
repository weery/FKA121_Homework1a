% plot the displacements
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('displacement.dat');

sizeD=size(data);

%plot
figure;


traj=randi([1, 256]);
traj_dat=zeros(1,10000);
for t=1:10000
        traj_dat(t)=norm([data(1,1+(traj-1)*3)-data(t,1+(traj-1)*3),data(1,2+(traj-1)*3)-data(t,2+(traj-1)*3),data(1,3+(traj-1)*3)-data(t,3+(traj-1)*3)]);
        traj_dat(t) = traj_dat(t) / 4.046;
end
plot(traj_dat)
% labels
xlabel('Time / [dim. unit]');
ylabel('Displacement / [dim. unit]');

% legend
legend('Atom 1','Atom 2','Atom 3');

% axis limits

