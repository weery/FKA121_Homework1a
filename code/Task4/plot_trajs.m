% Plot the displacements
% Made by Victor Nilsson and Simon Nilsson 2016-11-15

% load the data file
data = importdata('displacement.dat');

sizeD=size(data);

%plot
figure;

steplength = 1;
L = size(data);

traj=randi([1, (L(2)-1)/3], 3, 1);
traj_dat=zeros(3, L(1), 3);
for t=1:steplength:L(1)
    for i = 1:3
        traj_dat(:,t,i) = [...
            data(t,2+(traj(i)-1)*3);
            data(t,3+(traj(i)-1)*3);
            data(t,4+(traj(i)-1)*3)];
        % Divide by lattice parameter to obtain relative movement
        %traj_dat(:,t,i) = traj_dat(:,t,i) / 4.046;
    end
end

hold on
for i=1:3
   plot3(traj_dat(1,1:steplength:L(1),i), traj_dat(2,1:steplength:L(1),i), traj_dat(3,1:steplength:L(1),i)); 
end
hold off
%axis([0 17 0 17 0 17])

% labels
%xlabel('Time / [dim. unit]');
%ylabel('Displacement / [dim. unit]');

% legend
legend('Atom 1','Atom 2','Atom 3');

% axis limits

