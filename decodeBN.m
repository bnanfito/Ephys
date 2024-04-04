clear all
close all

dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
load(fullfile(dataFold,'dataSets/training/coolV1_MUdataSet.mat'))

for t = 1:4

    if t == 1
        tbl = v1bf{end}(v1bf{end}.goodUnit==1,:);
        tblLbl{t} = 'V1 Before';
    elseif t == 2
        tbl = v1af{end}(v1af{end}.goodUnit==1,:);
        tblLbl{t} = 'V1 After';
    elseif t == 3
        tbl = pssbf{end}(pssbf{end}.goodUnit==1,:);
        tblLbl{t} = 'PSS Before';
    elseif t == 4
        tbl = pssaf{end}(pssaf{end}.goodUnit==1,:);
        tblLbl{t} = 'PSS After';
    end

    obs = cat(3,tbl.fr.bc);
    obs = obs(:,1:end-1,:);
    id = repmat(1:size(obs,2),size(obs,1),1);
    Obs = reshape(obs,size(obs,1)*size(obs,2),size(obs,3));
    ID = reshape(id,size(id,1)*size(id,2),size(id,3));
    badObs = sum(isnan(Obs),2)>0;
%     Obs = Obs(~badObs,:);
%     ID = ID(~badObs,:);
    ids = unique(ID);
    colors = zeros(3,length(ids));
    stepsize = 2/(length(ids)-1);
    colors(1,:) = [(0:stepsize:1) fliplr(0:stepsize:1)];
    colors(3,:) = [flip(0:stepsize:1) (0:stepsize:1)];
%     Obs = zscore(Obs);
    
%     numTrees = 1000;
%     forest = TreeBagger(numTrees,Obs,ID,'OOBPrediction','on','Method','Classification');
%     % oobLabels = oobPredict(forest);
%     % for i = 1:length(oobLabels)
%     %     guess(i,1) = str2double(oobLabels{i});
%     % end
%     % ind = randsample(length(oobLabels),10);
%     % table(ID(ind),oobLabels(ind),...
%     %     VariableNames=["TrueLabel" "PredictedLabel"])
%     oobErr = oobError(forest);
%     oobAcc(t) = 1-oobErr(end);

    [coeff,score,latent,tsquare] = pca(Obs,'Economy',false);
    figure
    subplot(2,1,1);hold on
    for i = 1:length(ids)
        plot3(score(ID==ids(i),1),score(ID==ids(i),2),score(ID==ids(i),3),'o','Color',colors(:,i),'LineWidth',2)
        lbl{i} = [num2str(ids(i))]; 
    end
    legend(lbl)
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')

    subplot(2,1,2);hold on
    plot(latent,'-o','LineWidth',2)

    clear colors
end
