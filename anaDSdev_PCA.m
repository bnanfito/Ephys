clear all
close all

% dataFold = '/Volumes/Lab drive/Brandon/data/dataSets/DSdev';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/DSdev';
% dataFold = 'F:\Brandon\data\dataSets\DSdev';
dataFold = 'Y:\Brandon\data\dataSets\DSdev';
load(fullfile(dataFold,"DSdev_dataSet.mat"))

areas = {'V1','PSS'};
ageGroups = {[29 32],[33 36],[37 max(projectTbl.age)]};
% ageGroups = {[28 32],[29 33],[30 34],[31 35],[32 36],[33 37],[34 38],[35 39],[36 40],[37 41],[38 42],[39 43],[40 44],[41 300]};

nAG = length(ageGroups);
nAR = length(areas);
for ar = 1:nAR
for ag = 1:nAG
    ageLims = ageGroups{ag};
    areaIdx = strcmp(projectTbl.recSite,areas{ar});
    ageLimIdx = projectTbl.age>=ageLims(1) & projectTbl.age<=ageLims(2);
    dat{ar,ag} = vertcat(projectTbl(areaIdx & ageLimIdx,:).sumStats{:});

    dat{ar,ag} = dat{ar,ag}(dat{ar,ag}.goodUnit,:);

    [~,sortIdx] = sort(dat{ar,ag}.oriPref);
    dat{ar,ag} = dat{ar,ag}(sortIdx,:);
    nU(ar,ag) = height(dat{ar,ag});

    [distDat{ar,ag}] = anaPCA(dat{ar,ag});


%     [~] = popDecode(rTrial{ar,ag},cTrial{ar,ag});

end
end

%% PLOT

f = figure;
f.Position = [100 100 1000 600];
for ar = 1:nAR
    if strcmp(areas{ar},'V1')
        clr = 'b';
    elseif strcmp(areas{ar},'PSS')
        clr = 'r';
    end
    for ag = 1:nAG
    
        curDat = distDat{ar,ag};
        subplot(nAR,nAG,ag+(nAG*(ar-1)))
        
%         imagesc(curDat.rdm)
%         colorbar
%         if ag == 1
%             nYtick = length(curDat.cMean);
%             yticks(1:2:nYtick)
%             yticklabels(curDat.cMean(1:2:nYtick))
%         else
%             yticks([])
%         end
%         if ar == nAR
%             nXtick = length(curDat.cMean);
%             xticks(1:2:nXtick)
%             xticklabels(curDat.cMean(1:2:nXtick))
%         else
%             xticks([])
%         end
        
%         hold on
%         score = vertcat(curDat.score,curDat.score(1,:));
%         plot3(score(:,1),score(:,2),score(:,3),'k','LineWidth',1.5)
%         np = size(curDat.score,1);
%         dirClrs = hsv(np);
%         oriClrs = repmat(hsv(np/2),2,1);
%         for i = 1:np
%             plot3(score(i,1),score(i,2),score(i,3),'.','Color',dirClrs(i,:),'MarkerSize',30)
%         end
%         xlabel('PC1')
%         ylabel('PC2')
%         zlabel('PC3')
%         view(2)
%         oriID = mod(curDat.cMean,180)';
%         oris = unique(oriID);
%         for o = 1:length(oris)
%             plot3(curDat.score(oriID==oris(o),1),curDat.score(oriID==oris(o),2),curDat.score(oriID==oris(o),3),'--','Color',[0.6 0.6 0.6],'LineWidth',1.5)
%         end

        hold on
        dispIdx = curDat.cMean>0 & curDat.cMean<=180;
        angDisp = curDat.cMean(dispIdx);
        distF = curDat.distF(dispIdx);
        distNull = curDat.distNull(:,dispIdx);
        plot(angDisp,distF,'k-o','LineWidth',2)
        sem = std(curDat.rdmShift)/sqrt(size(curDat.rdmShift,1)); sem = sem(dispIdx);
        v = var(curDat.rdmShift); v = v(dispIdx);
        patch([angDisp fliplr(angDisp)],[distF-sem fliplr(distF+sem)],'k','EdgeColor','none','FaceAlpha',0.2)
        plot(angDisp,mean(distNull),'r-o')
        sem = std(distNull)/sqrt(size(distNull,1));
        v = var(distNull);
        sig2 = std(distNull,'omitnan')*2;
        for i = 1:size(distNull,2)
            P99(i) = prctile(distNull(:,i),99,'Method','exact');
            P95(i) = prctile(distNull(:,i),95,'Method','exact');
            P05(i) = prctile(distNull(:,i),5,'Method','exact');
            P01(i) = prctile(distNull(:,i),1,'Method','exact');
        end
        patch([angDisp fliplr(angDisp)],[mean(distNull)-sem fliplr(mean(distNull)+sem)],'r','EdgeColor','none','FaceAlpha',0.2)
        patch([angDisp fliplr(angDisp)],[P05 fliplr(P95)],'r','FaceColor','none','EdgeColor','r','LineStyle','--')
        patch([angDisp fliplr(angDisp)],[P01 fliplr(P99)],'r','FaceColor','none','EdgeColor','r','LineStyle',':')
        sigHiX = angDisp(distF>P95);
        sigHiY = distF(distF>P95);
        sigLoX = angDisp(distF<P05);
        sigLoY = distF(distF<P05);
        text(sigHiX,sigHiY+(sigHiY*0.1),'*')
        text(sigLoX,sigLoY-(sigLoY*0.1),'*')
        xticks([min(angDisp) 90 max(angDisp)])
        xlim([min(angDisp) max(angDisp)])
        if ar == nAR
            xlabel('angular disparity (+/- deg)')
        end
        if ag == 1
            ylabel('Euclidean distance')
        end

        title([areas{ar} '; P' num2str(ageGroups{ag}(1)) '-' num2str(ageGroups{ag}(2))])
    end
end

