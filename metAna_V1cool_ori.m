%metAna_V1cool_ori
clear all
close all

load('/Volumes/Lab drive/Brandon/data/dataSets/cooling/V1cool_MU_ori/V1cool_MU_ori_projectTbl.mat');

pssIdx = strcmp(projectTbl.recSite,'PSS');
v1Idx = strcmp(projectTbl.recSite,'V1');

data = projectTbl(pssIdx,:);
animals = unique(data.experimentId);
for a = 1:length(animals)
    ages(a) = unique(data.age(strcmp(data.experimentId,animals{a})));
end

figure;hold on
bins = 0:1:100;
histogram(ages,bins)
