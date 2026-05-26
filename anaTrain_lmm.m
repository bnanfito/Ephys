%anaTrain_lmm
clear all
close all

load('/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/training/Train_V1Cool/MU/Train_V1Cool_MUdataSet.mat')

d1 = data.v1bf; d1 = d1(screenUnits(d1,'MU'),:); 
d2 = data.v1af; d2 = d2(screenUnits(d2,'MU'),:);
d3 = data.pssbf; d3 = d3(screenUnits(d3,'MU'),:);
d4 = data.pssaf; d4 = d4(screenUnits(d4,'MU'),:);
D = vertcat(d1,d2,d3,d4);
priorFlag = [zeros(height(d1),1);ones(height(d2),1);zeros(height(d3),1);ones(height(d4),1)];
recSite = [repmat({'V1'},height(d1),1);repmat({'V1'},height(d2),1);repmat({'PSS'},height(d3),1);repmat({'PSS'},height(d4),1)];
trainType = repmat({'V1Cool'},height(D),1);

aniId = [];
for i = 1:height(D)
    aniId{i} = D.exptName{i}(1:5);
    age(i) = unique(projectTbl.age( strcmp(projectTbl.experimentId,aniId{i}) ));
end
D.animal = categorical(aniId');
D.priorFlag = priorFlag;
D.recSite = recSite;
D.trainType = trainType;
D.age = age';

exclude = {'febk7','febm6'};
D = D(~ismember(D.animal,exclude),:);

muTraining = table(D.animal,D.recSite,D.trainType,D.priorFlag,D.ldr,D.lor,D.age,'VariableNames',{'animal','recSite','trainType','priorFlag','Ldir','Lori','age'});


areas = unique(D.recSite);
nAR = length(areas);
figure;
for ar = 1:nAR
    subplot(1,2,ar);hold on
    if strcmp(areas{ar},'PSS')
        colororder([[0.1 0.5 1];[0 0 1]])
    elseif strcmp(areas{ar},'V1')
        colororder([[1 0.5 0.1];[1 0 0]])
    end
    arIdx = strcmp(D.recSite,areas{ar});
    Dar = D(arIdx,:);
    b = boxchart(Dar.animal,Dar.ldr,'GroupByColor',Dar.priorFlag,'Notch','on','MarkerStyle','.','MarkerSize',5);
    title(areas{ar})
    xlabel('animal')
    ylabel('Ldir')
    ylim([0 0.7])
    legend
end

figure;
for ar = 1:nAR
    subplot(1,2,ar);hold on
    if strcmp(areas{ar},'PSS')
        clr = 'r';
    elseif strcmp(areas{ar},'V1')
        clr = 'b';
    end
    arIdx = strcmp(D.recSite,areas{ar});
    Dar = D(arIdx,:);
    b = boxchart(Dar.priorFlag,Dar.ldr,'Notch','on','MarkerStyle','none','MarkerSize',5,'Color');
    b.JitterOutliers = 'on';
    b.MarkerStyle = '.';
    swarmchart(Dar.priorFlag,Dar.ldr,'XJitterWidth',0.5)
    title(areas{ar})
    xlabel('training')
    ylabel('Ldir')
    ylim([0 0.7])
    legend
end


% clear aniId
% 
% bf = data.pssbf; bf = bf(screenUnits(bf,'MU'),:);
% af = data.pssaf; af = af(screenUnits(af,'MU'),:);
% D = [bf;af];
% for i = 1:height(D)
%     aniId{i} = D.exptName{i}(1:5);
% end
% aniId = aniId';
% tPhase = [zeros(height(bf),1);ones(height(af),1)];
% 
% figure;
% boxchart(tPhase,D.ldr,'GroupByColor',aniId)

