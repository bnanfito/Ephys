%anaMerge
clear all
close all

dataFold = 'Y:\Brandon\data';

%% list animals

% load(fullfile(dataFold,'dataSets','cooling','V1cool_MU_ori','V1cool_MU_ori_projectTbl.mat'))
% for a = 1:length(animals)
%     cd(fullfile(dataFold,'Ephys',animals{a}))
%     folders = dir;
%     fileIdx = find(contains({folders.name},'MMM'));
%     if isempty(fileIdx)
%         continue
%     end
%     mergeName{a,1} = folders(fileIdx).name;
%     mergeId = mergeName{a}(12:end);
%     load(fullfile(dataFold,'Ephys',animals{a},mergeName{a},[mergeName{a} '_id.mat']))
%     probeId = find(strcmp({id.probes.area},'PSS'));
%     [sumStats{a}] = plotMerge(animals{a}, mergeId, probeId, dataFold, 0);
% 
% end

load(fullfile(dataFold,'dataSets','cooling','V1cool_MU_ori','anaMerge_dataSet.mat'))

LDR = [];
RP = [];
AGE = [];
for a = 1:length(animals)
    if isempty(sumStats{a})
        continue
    end
    dat = sumStats{a};
    nU = height(dat{1});

    cntrlLDR = dat{1}.ldr;
    coolLDR = dat{2}.ldr;
    LDR = vertcat(LDR,[cntrlLDR coolLDR]);

    cntrlrPref = dat{1}.rPref;
    coolrPref = dat{2}.rPref;
    RP = vertcat(RP,[cntrlrPref coolrPref]);

    AGE = vertcat(AGE,repmat(ages(a),nU,1));

end

figure; hold on
LDRdiff = (LDR(:,2)-LDR(:,1))./(LDR(:,2)+LDR(:,1));
plot(AGE,LDRdiff,'k.')
uAges = unique(ages);
for a = 1:length(uAges)
    x(a) = mean(LDRdiff(AGE==uAges(a)),'omitnan');
end
plot(uAges,x,'ko')
xlabel('age')
ylabel('delta Ldir')
clear x

figure;hold on
ageGroups = {[28 33],[34 40],[41 300]};
clrs = {[1 0 0], [0.5 0 0], [0 0 0]};
plot([0 1],[0 1],'k')
for ag = 1:length(ageGroups)
    idx = AGE>=ageGroups{ag}(1) & AGE<=ageGroups{ag}(2);
    plot(LDR(idx,1),LDR(idx,2),'.','Color',clrs{ag},'MarkerSize',10)
end
xlabel('control LDR')
ylabel('cooled LDR')



figure;hold on
rPrefDiff = (RP(:,2) - RP(:,1))./(RP(:,2) + RP(:,1));
plot(AGE,rPrefDiff,'k.')
for a = 1:length(uAges)
    x(a) = mean(rPrefDiff(AGE==uAges(a)),'omitnan');
end
plot(uAges,x,'ko')
xlabel('age')
ylabel('SI')











































































% % % clear all
% % % % close all
% % % 
% % % if ispc
% % % %     dataFold = 'D:\data'; 
% % % %     dataFold = 'C:\Users\brand\Documents\data';
% % % %     dataFold = 'F:\Brandon\data';
% % % dataFold = 'F:\Brandon\VSS2024\data';
% % % elseif ismac
% % % %     dataFold = '/Volumes/Lab drive/Brandon/data';
% % %     dataFold = '/Volumes/Lab drive/Brandon/VSS2024/data';
% % % %     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
% % % end
% % % physDir = fullfile(dataFold,'Ephys');
% % % dsDir = fullfile(dataFold,'dataSets');
% % % 
% % % anaMode = 'SU';
% % % for i = 1:3
% % %     figure; hold on
% % %     if i == 1
% % %         load(fullfile(dataFold,'dataSets',['P33to34Merge' anaMode 'Dat.mat']))
% % %         clr = [0.9290 0.6940 0.1250]; %yellow
% % %         ttl = 'P33-34';
% % %     elseif i == 2
% % %         load(fullfile(dataFold,'dataSets',['P36to37Merge' anaMode 'Dat.mat']))
% % %         clr = [0.8500 0.3250 0.0980]; %orange
% % %         ttl = 'P36-37';
% % %     elseif i == 3
% % %         load(fullfile(dataFold,'dataSets',['adultMerge' anaMode 'Dat.mat']))
% % %         clr = [0.6350 0.0780 0.1840]; %red
% % %         ttl = 'adult';
% % %     end
% % %     fileID = vertcat(dat.tbl.fileID{:});
% % %     condID = fileID(:,2);
% % %     idxA = dat.tbl.goodUnit == 1 & condID == 1;
% % %     idxB = dat.tbl.goodUnit == 1 & condID == 2;
% % % 
% % %     A = dat.tbl(idxA,:).DSI;
% % %     B = dat.tbl(idxB,:).DSI;
% % % 
% % %     A(A>1) = 1;
% % %     B(B>1) = 1;
% % % 
% % %     cdfA = cdfplot(A);
% % %     cdfA.Color = 'k';
% % %     
% % %     cdfB = cdfplot(B);
% % %     cdfB.Color = 'c';
% % % 
% % %     xlim([0 1])
% % %     title(ttl)
% % %     clear dat
% % % end