clear all
close all

dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
load(fullfile(dataFold,'dataSets','training','cntrlAL_MUdataSet.mat'))

dsFact = [1 2 3 4 5];
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};

figure
for t = 1:4

    if t == 1
        tbl = v1bf{end};
        tblLbl{t} = 'V1 Before';
    elseif t == 2
        tbl = v1af{end};
        tblLbl{t} = 'V1 After';
    elseif t == 3
        tbl = pssbf{end};
        tblLbl{t} = 'PSS Before';
    elseif t == 4
        tbl = pssaf{end};
        tblLbl{t} = 'PSS After';
    end
tbl = tbl(tbl.goodUnit==1,:); %exclude bad units

for ds = 1:length(dsFact)

    subTbl{t,ds} = [];
    expts = unique(tbl.exptName);
    for e = 1:length(expts)
        exptInd = find(strcmp(tbl.exptName,expts{e}));
        dsInd = downsample(exptInd,dsFact(ds));
        nU{t,ds,e} = length(dsInd);
        subTbl{t,ds} = vertcat(subTbl{t,ds},tbl(dsInd,:));
    end
    
    lbl{ds} = ['ds = ' num2str(ds) '; n = ' num2str(height(subTbl{t,ds}))];
    
    for e = 1:length(expts)
        exptInd = find(strcmp(subTbl{t,ds}.exptName,expts{e}));
        D{t,ds,e} = pdist([vertcat(subTbl{t,ds}(exptInd,:).xPos) vertcat(subTbl{t,ds}(exptInd,:).zPos)]);
        D{t,ds,e} = squareform(D{t,ds,e});
        for c = 1:size(D{t,ds,e},2)
            [neighD{t,ds,e}(:,c),neighN{t,ds,e}(:,c)] = sort(D{t,ds,e}(:,c));
            deltaDSI{t,ds,e}(:,c) = subTbl{t,ds}(exptInd,:).dsi(c)-subTbl{t,ds}(exptInd,:).dsi(neighN{t,ds,e}(:,c));
            absDeltaDSI{t,ds,e}(:,c) = abs(subTbl{t,ds}(exptInd,:).dsi(c)-subTbl{t,ds}(exptInd,:).dsi(neighN{t,ds,e}(:,c)));
        end
        if size(D{t,ds,e},2)<max([nU{t,ds,:}])
            deltaDSI{t,ds,e}(size(D{t,ds,e},2):max([nU{t,ds,:}]),:)=nan;
            absDeltaDSI{t,ds,e}(size(D{t,ds,e},2):max([nU{t,ds,:}]),:)=nan;
        end
        nearestD{t,ds,e} = neighD{t,ds,e}(2,:);

    end

    subplot(4,4,1+(4*(t-1)));hold on
    p = cdfplot([nearestD{t,ds,:}]);
    p.LineWidth = 2;
    xlabel('nearest neighbor distance (um)')
    ylabel([tblLbl{t} ' percentile'])
    if ds == length(dsFact)
        legend(lbl)
    end
    title([])

    subplot(4,4,2+(4*(t-1)));hold on  
    p = cdfplot(subTbl{t,ds}.dsi);
    p.LineWidth = 2;
    xlabel('dsi')
    ylabel('percentile')
    title([])

    if ds ==1
        subplot(4,4,3+(4*(t-1)));hold on
        y = mean([deltaDSI{t,ds,:}],2,'omitnan');
        plot(y,'-o','LineWidth',2)
        plot(repmat(1:length(y),2,1),y'+([1;-1]*std([deltaDSI{t,ds,:}],[],2,'omitnan')'),'Color',colors{ds},'LineWidth',2)
        xlabel('sorted neighbors')
        ylabel('delta DSI')
        ylim([-1 1])
        title([])
    end

    if ds == 1
        subplot(4,4,4+(4*(t-1)));hold on
        y = mean([absDeltaDSI{t,ds,:}],2,'omitnan');
        plot(y,'-o','LineWidth',2)
        plot(repmat(1:length(y),2,1),y'+([1;-1]*std([absDeltaDSI{t,ds,:}],[],2,'omitnan')'),'Color',colors{ds},'LineWidth',2)
        xlabel('sorted neighbors')
        ylabel('abs(delta DSI)')
        ylim([0 1])
        title([])
    end

end

end

