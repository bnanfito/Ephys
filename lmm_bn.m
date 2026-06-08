% linear mixed effects model - Brandon Nanfito
clear
% close all

% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
dataFold = 'Y:\Brandon\data';
proj = 'Train_V1Cool';
anaMode = 'MU';
area = 'PSS';
load(fullfile(dataFold,'dataSets','training',proj,anaMode,[proj '_' anaMode 'dataSet.mat']))

%% Organize data

d1 = data.v1bf; d1 = d1(screenUnits(d1,'MU'),:);
d2 = data.v1af; d2 = d2(screenUnits(d2,'MU'),:);
d3 = data.pssbf; d3 = d3(screenUnits(d3,'MU'),:);
d4 = data.pssaf; d4 = d4(screenUnits(d4,'MU'),:);
D_tmp = vertcat(d1,d2,d3,d4);
manipFlag = [zeros(height(d1),1);ones(height(d2),1);...
             zeros(height(d3),1);ones(height(d4),1)];
aniId{height(D_tmp),1} = [];
age = nan(height(D_tmp),1);
for i = 1:height(D_tmp)
    aniId{i} = D_tmp.exptName{i}(1:5);
    age(i) = unique(projectTbl.age( strcmp(projectTbl.experimentId,aniId{i}) ));
end
D_tmp.animalId = aniId;
D_tmp.age = age;
D_tmp.manipulation = manipFlag;
exclude = {'febk7','febm6'};
D_tmp = D_tmp(~ismember(D_tmp.animalId,exclude),:);
D_tmp = D_tmp(strcmp(D_tmp.area,area),:);
D = table(D_tmp.animalId, D_tmp.area, D_tmp.age, D_tmp.manipulation, D_tmp.ldr,...
          D_tmp.lor,'VariableNames',{'animalId','area','age','manipulation','ldr','lor'});



