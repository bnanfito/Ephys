%anaV1cool_contrast_merge
clear all
close all

proj = 'V1cool_contrast';
% dataFold = '/Volumes/NielsenHome2/Brandon/data';
dataFold = 'Y:\Brandon\data';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
anaMode = 'SU';

%% Build Project Table

projTbl = getProjectFiles(proj,1,'age','recSite','penNr','priorMFlag','priorDescr',...
                                       'duringMFlag','manipDescr','manipDetail',...
                                       'looperNameCond1','looperNameCond2',...
                                       'looperNameCond3','looperNameCond4',...
                                       'looperNameCond5');
animals = unique(projTbl.experimentId);
sumStats = [];
for a = 1:length(animals)
    aniIdx = strcmp(projTbl.experimentId,animals{a});
    ages(a) = unique(projTbl.age(aniIdx));
    cd(fullfile(dataFold,'Ephys',animals{a}))
    folders = dir;
    fileIdx = find(contains({folders.name},'MMM'));
    if isempty(fileIdx)
        sumStats = vertcat(sumStats,{[],[],[]});
        continue
    end

    out = cell(1,3);
    for fId = 1:length(fileIdx)
        mergeName{a,fId} = folders(fileIdx(fId)).name;
        mergeId = mergeName{a,fId}(12:end);
        load(fullfile(dataFold,'Ephys',animals{a},mergeName{a,fId},[mergeName{a,fId} '_id.mat']))
        probeId = find(strcmp({id.probes.area},'PSS'));
        disp(['generating sumStats for ' mergeName{a,fId}])
    
%         splitIntan(fullfile(dataFold,'Ephys'),animals{a},mergeId,probeId,'BRN')
        %Load merge info
        load(fullfile(dataFold,'Ephys',animals{a},mergeName{a,fId},[mergeName{a,fId} '_mergeInfo.mat']))
        nFiles = length(mergeInfo.files);
        for f = 1:nFiles
            exptName{f,1} = [animals{a} '_' mergeInfo.files{f}];
            load(fullfile(dataFold,'Ephys',animals{a},exptName{f},[exptName{f} '_trialInfo.mat']))
            isCon = length(trialInfo.dom)==2 & sum(strcmp(trialInfo.dom,'ori'))==1 & sum(strcmp(trialInfo.dom,'contrast'))==1 ;
            if isCon
                out_tmp = anaCon(animals{a},exptName{f}(8:10),exptName{f}(12:14),probeId,anaMode,dataFold,0,f);
                out{f} = vertcat(out{f},out_tmp);
            end
        end
    end

%     [out] = plotMerge(animals{a}, mergeId, probeId, dataFold, 0);

    sumStats = vertcat(sumStats,out);
end

% load(fullfile(dataFold,'dataSets','cooling',proj,'matchedSU',[proj '_matchedSUdataSet.mat']))

%% Organize Data

% compile sumStat tables by phase of experiment (1:cntrl, 2:cool, 3:post)
dat{1} = vertcat(sumStats{:,1});
dat{2} = vertcat(sumStats{:,2});
dat{3} = vertcat(sumStats{:,3});

% determine which units pass inclusion criteria for each phase of the expt
goodId(:,1) = screenUnits(dat{1},'SU');
goodId(:,2) = screenUnits(dat{2},'SU');
goodId(:,3) = screenUnits(dat{3},'SU');

keepIdx = goodId(:,1)&goodId(:,2);
dat{1} = dat{1}(keepIdx,:);
dat{2} = dat{2}(keepIdx,:);
dat{3} = dat{3}(keepIdx,:);

% for each phase, extract basic information like age, tuning metrics, etc
% and store that info in vectors
nU = height(dat{1});
uAge = nan(nU,1);
late = nan(nU,3);
rPref = nan(nU,3);

for u = 1:nU
    uAge(u) = ages(strcmp(animals,dat{1}.exptName{u}(1:5)));
end
for e = 1:3
%     ldr(:,e) = dat{e}.ldr;
    late(:,e) = dat{e}.latency;
%     rPref(:,e) = dat{e}.rPref;
    for u = 1:nU

        rMat = dat{e}.response{u};
        rMat(rMat<0) = 0;
        conds = dat{e}.condition{u}(strcmp(dat{e}.paramKey{u},'ori'),:);


        %tuning curve


        %PSTH
        binSize = 0.02;
        bins = -1:binSize:2;
        spkTimes = dat{e}.spkTimes{u}(1,:);
        nTrials = max(dat{e}.fr(u).trialNum,[],'all');
        psth(u,:,e) = histcounts(spkTimes,bins)/(nTrials*binSize);

    end
end



%% Plot


