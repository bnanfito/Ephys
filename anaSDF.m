clear all
close all

anaMode = 'MU';
% load(['/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/training/Train_V1Cool/' anaMode '/ranksum & rPref above 2/Train_V1Cool_' anaMode 'dataSet.mat'])

tbl = data.v1bf;
tbl = tbl(tbl.goodUnit,:);

[~,lateSortIdx] = sort(tbl.latency);
for u = 1:height(tbl)
    x = -1:0.01:2;
    spks = tbl.spkTimes{u}(1,:);
    g = normpdf(repmat(x,length(spks),1)',spks,repmat(0.05,1,length(spks)));
    sdf(u,:) = mean(g,2,'omitnan');
    sdf(u,:) = sdf(u,:)/max(sdf(u,:));
end

figure; hold on
imagesc(sdf(lateSortIdx,:))




