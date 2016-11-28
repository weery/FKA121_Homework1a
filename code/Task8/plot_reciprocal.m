% Plot the static structure factor using data from
% reciprocal lattice vectors
clear all

data = importdata('sq.dat');

% plot
figure;

% Remove the delta-distribution at distance zero
z=find(data(:,3)==0);
data(z,:)=[];

% Sort data after distance to centre
[sortedData,sortedIndex]=sort(data(:,3));
sortedData(:,2)=data(sortedIndex,2)+data(sortedIndex,1);

nBins=150;

maxBin = max(data(1,:));
minBin = min(data(1,:));
deltaBin = (maxBin-minBin)/nBins;

% Average in all radial bins
totBin=zeros(2,nBins);
idx=1;
for i=1:nBins
    totBin(1,i)=minBin+deltaBin*i;
    nCurrentBin=0;
    currentQ=sortedData(idx,1);

    while currentQ < totBin(1,i)
        nCurrentBin = nCurrentBin+1;
        totBin(2,i)=totBin(2,i)+sortedData(idx,2);
        idx=idx+1;
        currentQ=sortedData(idx,1);
    end
    totBin(2,i)= totBin(2,i)/nCurrentBin;
end

plot(totBin(1,:),totBin(2,:),'-')

xlim([0, totBin(1,end)])

title('Simulated Static structure factor','interpreter','latex')
xlabel('q','interpreter','latex')
ylabel('S(q)','interpreter','latex')