%plot LMM data based on saved table

%load data
load('lmm.mat'); %LMM

%% plot Ldir - separate plots for training type and area

trType={'V1Cool'};
recSite={'V1';'PSS'};

offsetGroup=2; %offset between pre and post
randWidth=0.5; %random shuffle width

for t=1:length(trType)
    for r=1:length(recSite)
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

%% same plots, but for Lori

trType={'V1Cool'};
recSite={'V1';'PSS'};

offsetGroup=2; %offset between pre and post
randWidth=0.5; %random shuffle width

for t=1:length(trType)
    for r=1:length(recSite)
        figure
        trainingLMM = trainLMMLori;
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
                dataPoints(count:count+length(idx)-1,2)=MUdata.Lori(idx);
                dataPoints(count:count+length(idx)-1,3)=a;
                count=count+length(idx);

            end
        end
        scatter(dataPoints(:,1),dataPoints(:,2),[],dataPoints(:,3))

        %add group means - the values for these are constant across the
        %table, so can just use first row
        hold on
        preMean=LMMdata.interLori(1);
        postMean=preMean+LMMdata.slopeLori(1);
        plot([1-randWidth Nanim+randWidth],[preMean preMean],'k-');
        plot([1-randWidth Nanim+randWidth]+Nanim+offsetGroup,[postMean postMean],'k-');

        %add animal means
        for a=1:Nanim
            plot([-randWidth +randWidth]+a,[LMMdata.animPreMeanLori(a) LMMdata.animPreMeanLori(a)],'r-');
            plot([-randWidth +randWidth]+a+Nanim+offsetGroup,[LMMdata.animPostMeanLori(a) LMMdata.animPostMeanLori(a)],'r-');
        end

        %label axes and plot
        xt1=[1:Nanim];
        xt2=[1:Nanim]+Nanim+offsetGroup;
        set(gca,'XTick',[xt1 xt2])
        set(gca,'Box','on')
        ylabel('Lori')
        title([trType{t} ', ' recSite{r}])
        pbaspect([0.9 1 1])
    end

end




