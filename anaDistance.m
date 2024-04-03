clear all
close all

dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
load(fullfile(dataFold,'dataSets','training','cntrlAL_MUdataSet.mat'))

figure
for t = 1:4

    if t == 1
        tbl = v1bf{end};
    elseif t == 2
        tbl = v1af{end};
    elseif t == 3
        tbl = pssbf{end};
    elseif t == 4
        tbl = pssaf{end};
    end
tbl = tbl(tbl.goodUnit==1,:); %exclude bad units

dsFact = [1 2 3 4 5];
for ds = 1:length(dsFact)

    subTbl{ds} = [];
    expts = unique(tbl.exptName);
    for e = 1:length(expts)
        exptInd = find(strcmp(tbl.exptName,expts{e}));
        dsInd = downsample(exptInd,dsFact(ds));
        subTbl{ds} = vertcat(subTbl{ds},tbl(dsInd,:));
    end
    
    lbl{ds} = ['ds = ' num2str(ds) '; n = ' num2str(height(subTbl{ds}))];
    
    for e = 1:length(expts)
        exptInd = find(strcmp(subTbl{ds}.exptName,expts{e}));
        D{ds,e} = pdist([vertcat(subTbl{ds}(exptInd,:).xPos) vertcat(subTbl{ds}(exptInd,:).zPos)]);
        D{ds,e} = squareform(D{ds,e});
        D{ds,e}(D{ds,e}==0)=nan;
        nearestD{ds,e} = min(D{ds,e});
    end

    subplot(4,4,1+(4*(t-1)));hold on
    cdfplot([nearestD{ds,:}]);
    xlabel('nearest neighbor distance (um)')
    ylabel('percentile')
    if ds == length(dsFact)
        legend(lbl)
    end
    title([])

    subplot(4,4,2+(4*(t-1)));hold on  
    cdfplot(subTbl{ds}.dsi)
    xlabel('dsi')
    title([])

    subplot(4,4,3+(4*(t-1)));hold on
    cdfplot(1-subTbl{ds}.dcv)
    xlabel('1-dcv')
    title([])

    subplot(4,4,4+(4*(t-1)));hold on
    cdfplot(subTbl{ds}.rPref)
    xlabel('rPref (ori)')
    title([])

end

end
