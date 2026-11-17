clear all
% close all

%% Load data

anaMode = 'SU';

% dataFold = '/Volumes/Lab drive/Brandon/data/dataSets/DSdev';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/DSdev';
% dataFold = '/Volumes/NielsenHome2/Brandon/data/dataSets/DSdev';
% dataFold = 'F:\Brandon\data\dataSets\DSdev';
dataFold = 'Y:\Brandon\data\dataSets\DSdev';
load(fullfile(dataFold,['DSdev_' anaMode 'dataSet.mat']))

areas = {'V1','PSS'};
ageGroups = {[29 32],[33 36],[37 max(projTbl.age)]};
% ageGroups = {[28 32],[29 33],[30 34],[31 35],[32 36],[33 37],[34 38],[35 39],[36 40],[37 41],[38 42],[39 43],[40 44],[41 300]};
nAG = length(ageGroups);
nAR = length(areas);

for ar = 1:nAR
for ag = 1:nAG
    
    ageLims = ageGroups{ag};
    areaIdx = strcmp(projTbl.recSite,areas{ar});
    ageLimIdx = projTbl.age>=ageLims(1) & projTbl.age<=ageLims(2);
    dat{ar,ag} = vertcat(projTbl(areaIdx & ageLimIdx,:).sumStats{:});

    %only take units that pass inclusion criteria (screenUnits)
    dat{ar,ag} = dat{ar,ag}(screenUnits(dat{ar,ag},anaMode),:);
    %sort units in data table by their preferred direction
    [~,sortIdx] = sort(dat{ar,ag}.oriPref);
    dat{ar,ag} = dat{ar,ag}(sortIdx,:);
    nU(ar,ag) = height(dat{ar,ag});

end
end


% oris = 0:90:270;
% oriColors = hsv(length(oris));
% nReps = 5;
% nTrials = nReps*length(oris);
% for ar = 1
% for ag = 3
%     D = dat{ar,ag};
%     for u = 1:nU(ar,ag)
% 
%         oriIdx = strcmp(D.paramKey{u},'ori');
%         cndIdx = ismember(D.cndKey{u}(:,oriIdx),oris);
%         sizeIdx = contains(D.paramKey{u},'size');
%         if length(D.paramKey{u})>1 && sum(sizeIdx)>0
%             cndIdx = cndIdx & D.cndKey{u}(:,sizeIdx)>100;
%         end
%         trials = D.fr(u).trialNum(1:nReps,cndIdx);
%         trials = trials(:)';
%         oriLbl = D.fr(u).trialCond(1:nReps,cndIdx);
%         oriLbl = D.cndKey{u}(oriLbl(:),oriIdx)';
%         repLbl = repmat(1:nReps,1,length(oris));
%         
%         for t = 1:length(trials)
%                 spks = D.spkTimes{u}(1,D.spkTimes{u}(2,:)==trials(t)); %find the spikes times for that trial...
%                 [sdfTrial{ar,ag}(:,t,u),time] = compute_sdf2(spks,0.001,[-1 1]); %and make a spike density function from those spike times
%         end
% 
%         oriLbl = repmat(oriLbl,length(time),1);
%         repLbl = repmat(repLbl,length(time),1);
%         timeLbl = repmat(time',1,nTrials);
% 
%     end
%     R = reshape(sdfTrial{ar,ag},nTrials*length(time),nU(ar,ag));
%     R = R(:,~sum(R)==0); %remove units that do not spike for these oris
%     R = R./max(R); %normalize the response rate for each unit (between 0 & 1)
%     oriLbl = oriLbl(:);
%     repLbl = repLbl(:);
%     timeLbl = timeLbl(:);
% 
%     [c,s] = pca(R);
% 
%     figure;hold on
%     for o = 1:length(oris)
%         oIdx = oriLbl==oris(o);
%         for r = 1:nReps
%             rIdx = oIdx&repLbl==r;
% 
%             plot3(s(rIdx,4),s(rIdx,2),s(rIdx,3),'Color',oriColors(o,:))
% 
%         end
%     end
% 
% end
% end




oris = 0:90:270;
oriColors = hsv(length(oris));
for ar = 1
for ag = 2
    for u = 1:nU(ar,ag)
        for o = 1:length(oris)
            cndIdx = dat{ar,ag}.cndKey{u}(:,strcmp(dat{ar,ag}.paramKey{u},'ori'))==oris(o);
            if length(dat{ar,ag}.paramKey{u})>1 && sum(contains(dat{ar,ag}.paramKey{u},'size'))>0
                cndIdx = cndIdx & dat{ar,ag}.cndKey{u}(:,contains(dat{ar,ag}.paramKey{u},'size'))>100;
            end
            trials = dat{ar,ag}.fr(u).trialNum(:,cndIdx); %index trials of this condition (direction)
            for t = 1:length(trials) %for each trial of this condition...
                spks = dat{ar,ag}.spkTimes{u}(1,dat{ar,ag}.spkTimes{u}(2,:)==trials(t)); %find the spikes times for that trial...
                [sdfTrial(t,:),time] = compute_sdf2(spks,0.001,[-1 2]); %and make a spike density function from those spike times
            end
            oLbl_mean(:,o) = repmat(oris(o),length(time),1);
            sdfMean(:,o,u) = mean(sdfTrial,'omitnan');
            sdfSEM(:,o,u) = sem(sdfTrial);
            clear sdfTrial
        end
    end
    
    %create matrix R where each column is a unit and each row is an observation (point in time for a particular mean orientation sdf)
    R = reshape(sdfMean,length(oris)*length(time),nU(ar,ag)); 
    oLbl_mean = oLbl_mean(:);
    R = R(:,~sum(R)==0); %remove units that do not spike for these oris
    R = R./max(R); %normalize the response rate for each unit (between 0 & 1)
    [c,s,~,~,e] = pca(R); %project data onto PCs
    figure; hold on
    for o = 1:length(oris)
        oIdx = oLbl_mean==oris(o);
        plot3(s(oIdx,5),s(oIdx,4),s(oIdx,3),'Color',oriColors(o,:))
    end
end
end