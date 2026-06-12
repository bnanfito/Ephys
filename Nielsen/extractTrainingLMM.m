%mixed linear effects model for training data - save results in matlab file

%%

%anaTrain_lmm
clear all
close all

type='V1cool';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
dataFold = 'Y:\Brandon\data';
load(fullfile(dataFold,'dataSets','training',['Train_' type],'MU',['Train_' type '_MUdataSet.mat']))

d1 = data.v1bf; d1 = d1(screenUnits(d1,'MU'),:); 
d2 = data.v1af; d2 = d2(screenUnits(d2,'MU'),:);
d3 = data.pssbf; d3 = d3(screenUnits(d3,'MU'),:);
d4 = data.pssaf; d4 = d4(screenUnits(d4,'MU'),:);
D = vertcat(d1,d2,d3,d4);
priorFlag = [zeros(height(d1),1);ones(height(d2),1);zeros(height(d3),1);ones(height(d4),1)];
recSite = [repmat({'V1'},height(d1),1);repmat({'V1'},height(d2),1);repmat({'PSS'},height(d3),1);repmat({'PSS'},height(d4),1)];
trainType = repmat({type},height(D),1);

aniId = [];
for i = 1:height(D)
    aniId{i} = D.exptName{i}(1:5);
    age(i) = unique(projectTbl.age( strcmp(projectTbl.experimentId,aniId{i}) ));
end
D.animal = aniId';
D.priorFlag = priorFlag;
D.recSite = recSite;
D.trainType = trainType;
D.age = age';

exclude = {'febk7','febm6'};
D = D(~ismember(D.animal,exclude),:);

muTraining = table(D.animal,D.recSite,D.trainType,D.priorFlag,D.ldr,D.lor,D.age,'VariableNames',{'animal','recSite','trainType','priorFlag','Ldir','Lori','age'});

%%
% %load training data
% load('C:\Kristina\paper\Khamiss et al 24 ds dev\data\MUtraining_table_041425.mat')
% 
% %create one long table
% muTraining=[muBiDirTraining;muStaticTraining;muGrayTraining];
% muTraining.priorFlag=categorical(muTraining.priorFlag);
% muTraining.animal=categorical(muTraining.animal);

trType={type};
recSite={'V1';'PSS'};
curSite = 2;

% %% add age for each animal
% trAnim=unique(muTraining.animal);
% connSQL=database('ephysDatabase','','');
% for i=1:length(trAnim)
%     selquery=['SELECT * FROM tblanimal WHERE experimentId="' char(trAnim(i)) '"'];
%     data=select(connSQL,selquery);
% 
%     idx=find(muTraining.animal==trAnim(i));
%     muTraining.age(idx)=data.age;
% end
% close(connSQL);

%% loop through training type and area, compute training effects
trainLMMLdir=table;
% trainLMMLori=table;


for t=1:length(trType)
    for r=curSite

        idx=find(strcmp(muTraining.recSite,recSite{r}) & strcmp(muTraining.trainType,trType{t}));
        Nanim=length(unique(muTraining.animal(idx)));


        %Ldir data
        infoAnimal=table;

        Tdata=muTraining(idx,[1 4 5]);
        lme_2ML=fitlme(Tdata,'Ldir ~ 1+ priorFlag + (1+priorFlag|animal)','FitMethod','ML');

        %general slope
        lmmInter=lme_2ML.Coefficients.Estimate(1); %mean for all animals before training
        lmmSlope=lme_2ML.Coefficients.Estimate(2); %slope for all animals

        infoAnimal.slopeLdir(1:Nanim)=lmmSlope;
        infoAnimal.interLdir(1:Nanim)=lmmInter;

        %values per animal
        [LMMval,LMMinfo,LMMstats]=randomEffects(lme_2ML);

        %animal slopes
        idx2=strcmp([LMMinfo.Name],'priorFlag');
        infoAnimal.Names=LMMinfo.Level(idx2);
        infoAnimal.animSlopeLdir=LMMval(idx2); %animal-specific slope
        infoAnimal.animSELdir=LMMstats.SEPred(idx2); %SE for the slope

        %animal intercept
        idx2=strcmp([LMMinfo.Name],'(Intercept)');
        infoAnimal.animInterLdir=LMMval(idx2);  %animal specific pre-mean

        %general info
        infoAnimal.trainType(:)=trType(t);
        infoAnimal.recSite(:)=recSite(r);

        %age
        for i=1:Nanim
            idxA=find(strcmp(muTraining.animal,infoAnimal.Names{i}),1);
            infoAnimal.age(i)=muTraining.age(idxA);
        end

        trainLMMLdir=[trainLMMLdir;infoAnimal];
        
        % %Lori data
        % infoAnimal=table;
        % 
        % Tdata=muTraining(idx,[1 4 6]);
        % lme_2ML=fitlme(Tdata,'Lori ~ 1+ priorFlag + (1+priorFlag|animal)','FitMethod','ML');
        % 
        % lmmInter=lme_2ML.Coefficients.Estimate(1);
        % lmmSlope=lme_2ML.Coefficients.Estimate(2);
        % infoAnimal.slopeLori(1:Nanim)=lmmSlope;
        % infoAnimal.interLori(1:Nanim)=lmmInter;
        % 
        % [LMMval,LMMinfo,LMMstats]=randomEffects(lme_2ML);
        % 
        % %animal slope
        % idx2=strcmp([LMMinfo.Name],'priorFlag');
        % infoAnimal.Names=LMMinfo.Level(idx2);
        % infoAnimal.animSlopeLori=LMMval(idx2); 
        % infoAnimal.animSELdir=LMMstats.SEPred(idx2); %SE for the slope
        % 
        % %animal pre-mean
        % idx2=strcmp([LMMinfo.Name],'(Intercept)');
        % infoAnimal.animInterLori=LMMval(idx2);
        % 
        % %general info
        % infoAnimal.trainType(:)=trType(t);
        % infoAnimal.recSite(:)=recSite(r);
        % 
        % %age
        % for i=1:Nanim
        %     idxA=find(strcmp(muTraining.animal,infoAnimal.Names{i}),1);
        %     infoAnimal.age(i)=muTraining.age(idxA);
        % end
        % 
        % trainLMMLori=[trainLMMLori;infoAnimal];

    end
end


%% compute additional values per animal
%full slope per animal
trainLMMLdir.animFullSlopeLdir=trainLMMLdir.animSlopeLdir+trainLMMLdir.slopeLdir;
% trainLMMLori.animFullSlopeLori=trainLMMLori.animSlopeLori+trainLMMLori.slopeLori;

%pre mean per animal
trainLMMLdir.animPreMeanLdir=trainLMMLdir.interLdir+trainLMMLdir.animInterLdir;
% trainLMMLori.animPreMeanLori=trainLMMLori.interLori+trainLMMLori.animInterLori;

%post mean per animal
trainLMMLdir.animPostMeanLdir=trainLMMLdir.interLdir+trainLMMLdir.slopeLdir+...
    trainLMMLdir.animInterLdir+trainLMMLdir.animSlopeLdir;
% trainLMMLori.animPostMeanLori=trainLMMLori.interLori+trainLMMLori.slopeLori+...
%     trainLMMLori.animInterLori+trainLMMLori.animSlopeLori;



%% save

save(fullfile(dataFold,'lmm.mat'),'trainLMMLdir','muTraining');
% save(fullfile(dataFold,'lmm.mat'),'trainLMMLori','muTraining');

%% plot

% trType={'V1Cool'};
% recSite={'V1';'PSS'};
% curSite = 1;

offsetGroup=2; %offset between pre and post
randWidth=0.5; %random shuffle width

for t=1:length(trType)
    for r=curSite
        figure
        trainingLMM = trainLMMLdir;
        %find relevant LMM data
        idx=find(strcmp(trainingLMM.trainType,trType{t}) & strcmp(trainingLMM.recSite,recSite{r}));
        LMMdata=trainingLMM(idx,:);

        %find relevant raw data
        idx=find(strcmp(muTraining.trainType,trType{t}) & strcmp(muTraining.recSite,recSite{r}));
        MUdata=muTraining(idx,:);

        Uanim=unique(MUdata.animal);
        Nanim=length(Uanim);

        dataPoints=zeros(height(MUdata),3);

        %plot raw data
        count=1;
        for p=0:1
            for a=1:Nanim
                centerBin=a+p*(offsetGroup+Nanim);
                idx=find(strcmp(MUdata.animal,Uanim{a}) & MUdata.priorFlag==p);

                dataPoints(count:count+length(idx)-1,1)=centerBin+rand(length(idx),1)*randWidth-randWidth/2; %rand is between 0 and 1
                dataPoints(count:count+length(idx)-1,2)=MUdata.Ldir(idx);
                dataPoints(count:count+length(idx)-1,3)=a;
                count=count+length(idx);

            end
        end
        scatter(dataPoints(:,1),dataPoints(:,2),[],dataPoints(:,3))

        %add group means - the values for these are constant across the
        %table, so can just use first row
        hold on
        preMean=LMMdata.interLdir(1);
        postMean=preMean+LMMdata.slopeLdir(1);
        plot([1-randWidth Nanim+randWidth],[preMean preMean],'k-');
        plot([1-randWidth Nanim+randWidth]+Nanim+offsetGroup,[postMean postMean],'k-');

        %add animal means
        for a=1:Nanim
            plot([-randWidth +randWidth]+a,[LMMdata.animPreMeanLdir(a) LMMdata.animPreMeanLdir(a)],'r-');
            plot([-randWidth +randWidth]+a+Nanim+offsetGroup,[LMMdata.animPostMeanLdir(a) LMMdata.animPostMeanLdir(a)],'r-');
        end

        %label axes and plot
        xt1=[1:Nanim];
        xt2=[1:Nanim]+Nanim+offsetGroup;
        set(gca,'XTick',[xt1 xt2])
        set(gca,'Box','on')
        ylabel('Ldir')
        title([trType{t} ', ' recSite{r}])
        pbaspect([0.9 1 1])
    end

end

% %% same plots, but for Lori
% 
% % trType={'V1Cool'};
% % recSite={'V1';'PSS'};
% 
% offsetGroup=2; %offset between pre and post
% randWidth=0.5; %random shuffle width
% 
% for t=1:length(trType)
%     for r=curSite
%         figure
%         trainingLMM = trainLMMLori;
%         %find relevant LMM data
%         idx=find(strcmp(trainingLMM.trainType,trType{t}) & strcmp(trainingLMM.recSite,recSite{r}));
%         LMMdata=trainingLMM(idx,:);
% 
%         %find relevant raw data
%         idx=find(strcmp(muTraining.trainType,trType{t}) & strcmp(muTraining.recSite,recSite{r}));
%         MUdata=muTraining(idx,:);
% 
%         Uanim=unique(MUdata.animal);
%         Nanim=length(Uanim);
% 
%         dataPoints=zeros(height(MUdata),3);
% 
%         %plot raw data
%         count=1;
%         for p=0:1
%             for a=1:Nanim
%                 centerBin=a+p*(offsetGroup+Nanim);
%                 idx=find(strcmp(MUdata.animal,Uanim{a}) & MUdata.priorFlag==p);
% 
%                 dataPoints(count:count+length(idx)-1,1)=centerBin+rand(length(idx),1)*randWidth-randWidth/2; %rand is between 0 and 1
%                 dataPoints(count:count+length(idx)-1,2)=MUdata.Lori(idx);
%                 dataPoints(count:count+length(idx)-1,3)=a;
%                 count=count+length(idx);
% 
%             end
%         end
%         scatter(dataPoints(:,1),dataPoints(:,2),[],dataPoints(:,3))
% 
%         %add group means - the values for these are constant across the
%         %table, so can just use first row
%         hold on
%         preMean=LMMdata.interLori(1);
%         postMean=preMean+LMMdata.slopeLori(1);
%         plot([1-randWidth Nanim+randWidth],[preMean preMean],'k-');
%         plot([1-randWidth Nanim+randWidth]+Nanim+offsetGroup,[postMean postMean],'k-');
% 
%         %add animal means
%         for a=1:Nanim
%             plot([-randWidth +randWidth]+a,[LMMdata.animPreMeanLori(a) LMMdata.animPreMeanLori(a)],'r-');
%             plot([-randWidth +randWidth]+a+Nanim+offsetGroup,[LMMdata.animPostMeanLori(a) LMMdata.animPostMeanLori(a)],'r-');
%         end
% 
%         %label axes and plot
%         xt1=[1:Nanim];
%         xt2=[1:Nanim]+Nanim+offsetGroup;
%         set(gca,'XTick',[xt1 xt2])
%         set(gca,'Box','on')
%         ylabel('Lori')
%         title([trType{t} ', ' recSite{r}])
%         pbaspect([0.9 1 1])
%     end
% 
% end








