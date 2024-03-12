%process units
clear all
close all

visTest = 'anova'; alpha = 0.01;
svePlt = 0;
chkSum = 0;
sveSum = 0;

if ispc
    dataFold = 'D:\data'; 
%     dataFold = 'F:\Brandon\data';
elseif ismac
    dataFold = '/Volumes/Lab drive/Brandon/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');
figDir = fullfile(dataFold,'Figures');
sumDir = fullfile(dataFold,'SummaryStats');
anaMode = 'MU';

% V1 COOLED
animals = {'febh2','febh3','febj3'};

% % CONTROL
% animals = {'febj2','febj4'};

% % CONTROL (AUGUSTO)
% animals = {'FEAO4','FEAQ5','FEAS6','FEAT1','FEAN6'};

% % CONTROL GRAY SCREEN (AUGUSTO)
% animals = {'FEAQ2','FEAQ3','FEAQ4','FEAQ7'};

% % STATIC (AUGUSTO)
% animals = {'FEAS9','FEAT2','FEAU5','FEAU8','FEAU9','FEAV0','FEAV1'};

nAnimals = length(animals);
for a = 1:nAnimals
%% Define experiments and training parameters
if strcmp(animals{a},'febh2')
    before = {'_u000_001'};
    after = {'_u000_040','_u000_041','_u002_000'};
    trainAx = [90 270];
elseif strcmp(animals{a},'febh3')
    before = {'_u000_002'};
    after = {'_u000_043','_u001_000','_u002_002','_u003_000'};
    trainAx = [0 180];
elseif strcmp(animals{a},'febj3')
    before = {'_u000_001'};
    after = {'_u000_038'};
    trainAx = [0 180];
elseif strcmp(animals{a},'febj2')
    before = {'_u000_003'};
    after = {'_u000_023'};
    trainAx = [90 270];
elseif strcmp(animals{a},'febj4')
    before = {'_u000_006'};
    after = {'_u000_027'};
    trainAx = [90 270];
elseif strcmp(animals{a},'FEAO4')
    before = {'_u000_001'};
    after = {'_u002_000'};
    trainAx = [0 180];
elseif strcmp(animals{a},'FEAQ5')
    before = {'_u000_000'};
    after = {'_u002_000'};
    trainAx = [90 270];
elseif strcmp(animals{a},'FEAS6')
    before = {'_u000_003'};
    after = {'_u002_001'};
    trainAx = [60 240];
elseif strcmp(animals{a},'FEAT1')
    before = {'_u000_000'};
    after = {'_u002_000'};
    trainAx = [90 270];
elseif strcmp(animals{a},'FEAN6')
    before = {'_u000_000'};
    after = {'_u002_002'};
    trainAx = [150 330];
elseif strcmp(animals{a},'FEAQ2')
    before = {'_u000_001'};
    after = {'_u002_001'};
    trainAx = [];
elseif strcmp(animals{a},'FEAQ3')
    before = {'_u000_001'};
    after = {'_u001_001'};
    trainAx = [];
elseif strcmp(animals{a},'FEAQ4')
    before = {'_u000_001'};
    after = {'_u001_001'};
    trainAx = [];
elseif strcmp(animals{a},'FEAQ7')
    before = {'_u000_000'};
    after = {'_u002_001'};
    trainAx = [];
elseif strcmp(animals{a},'FEAS9')
    before = {'_u000_000'};
    after = {'_u002_001'};
    trainAx = [120];
elseif strcmp(animals{a},'FEAT2')
    before = {'_u000_000'};
    after = {'_u002_002'};
    trainAx = [60];
elseif strcmp(animals{a},'FEAU5')
    before = {'_u000_000'};
    after = {'_u002_002'};
    trainAx = [120];
elseif strcmp(animals{a},'FEAU8')
    before = {'_u000_002'};
    after = {'_u002_003'};
    trainAx = [0];
elseif strcmp(animals{a},'FEAU9')
    before = {'_u000_000'};
    after = {'_u002_003'};
    trainAx = [90];
elseif strcmp(animals{a},'FEAV0')
    before = {'_u000_000'};
    after = {'_u002_000'};
    trainAx = [60];
elseif strcmp(animals{a},'FEAV1')
    before = {'_u000_001'};
    after = {'_u002_002'};
    trainAx = [150];
end
if isempty(trainAx)
    nFig = 1;
else
    orthAx = sort(mod(trainAx+90,360));
    nFig = 3;
end

v1bf{a} = [];
v1af{a} = [];
pssbf{a} = [];
pssaf{a} = [];

%% Extract summary statistics 

for tr = 1:2 % tr = 1 = before; tr = 2 = after training
    if tr == 1
        exptList = before;
    elseif tr == 2
        exptList = after;
    end

    for e = 1:length(exptList)
        exptName = [animals{a} exptList{e}]; unit = exptName(8:10); expt = exptName(12:14);
        load(fullfile(physDir,animals{a},exptName,[exptName '_id.mat']))
        for p = 1:length(id.probes)
            sumFile = fullfile(sumDir,animals{a},exptName,[exptName '_p' num2str(p) '_sumStats' anaMode '.mat']);
            if isfile(sumFile) && chkSum == 1
                disp(['loading sumStats for ' exptName])
                load(sumFile,'sumStats')
            else
                disp(['generating sumStats for' exptName])
                sumStats = plotUnits(animals{a},unit,expt,p,anaMode,visTest,alpha,0,sveSum,0,dataFold);close all
            end
%             sumStats = convertAL(animals{a},unit,expt,p,dataFold);
%             sumStats = sumStats(sumStats.goodUnit == 1,:);
            sumStats.dsi(sumStats.dsi>1) = 1;
            sumStats.dcv(sumStats.dcv>1) = 1;
            sumStats.dcv(sumStats.dcv<0) = 0;
            sumStats.osi(sumStats.osi>1) = 1;

            if ~isempty(trainAx)
                oriDiff = abs(sumStats.cPref-trainAx);
                oriDiff(oriDiff>180) = 360-oriDiff(oriDiff>180);
                trainAxID = sum(oriDiff<=30,2)==1;
                sumStats.trainAx = trainAxID;
                
                oriDiff = abs(sumStats.cPref-orthAx);
                oriDiff(oriDiff>180) = 360-oriDiff(oriDiff>180);
                orthAxID = sum(oriDiff<=30,2)==1;
                sumStats.orthAx = orthAxID;
            end

            if strcmp(id.probes(p).area,'V1')
                if tr == 1
                    v1bf{a} = vertcat(v1bf{a},sumStats);
                elseif tr == 2
                    v1af{a} = vertcat(v1af{a},sumStats);
                end
            elseif strcmp(id.probes(p).area,'PSS')
                if tr == 1
                    pssbf{a} = vertcat(pssbf{a},sumStats);
                elseif tr == 2
                    pssaf{a} = vertcat(pssaf{a},sumStats);
                end
            end
        end
    end

end

end

animals{nAnimals+1} = 'all animals';
v1bf{nAnimals+1} = vertcat(v1bf{1:nAnimals});
v1af{nAnimals+1} = vertcat(v1af{1:nAnimals});
pssbf{nAnimals+1} = vertcat(pssbf{1:nAnimals});
pssaf{nAnimals+1} = vertcat(pssaf{1:nAnimals});

%% PLOT

metrics = {'rPref','dsi','dcv','osi'};
nMet = length(metrics);

for f = 1:nFig % all units; pref trained; pref orth
    figure
    if f == 1
        ttl{f} = 'all units';
    elseif f == 2
        ttl{f} = 'pref train';
    elseif f == 3
        ttl{f} = 'pref orth';
    end

    for a = 1:nAnimals+1
    for m = 1:nMet+1
        for t = 1:4 % each area/before and after
        
            if t == 1 % V1 before training
                if isempty(v1bf{a})
                    lbl{a,t,f} = '';
                    lblBit(t) = false;
                    continue
                end
                linStyl = '--';
                col = 'b';
                mrk = 'o';
                mrkSz = 5;
                if f == 1 % all units
                    tbl = v1bf{a};
                elseif f == 2 % pref train
                    tbl = v1bf{a}(v1bf{a}.trainAx==1,:);
                elseif f == 3 % pref orth
                    tbl = v1bf{a}(v1bf{a}.orthAx==1,:);
                end
                lbl{a,t,f} = 'V1 before; n=';
            elseif t == 2 % V1 after training
                if isempty(v1af{a})
                    lbl{a,t,f} = '';
                    lblBit(t) = false;
                    continue
                end
                linStyl = '-';
                col = 'b';
                mrk = '.';
                mrkSz = 20;
                if f == 1
                    tbl = v1af{a};
                elseif f == 2
                    tbl = v1af{a}(v1af{a}.trainAx==1,:);
                elseif f == 3
                    tbl = v1af{a}(v1af{a}.orthAx==1,:);
                end
                lbl{a,t,f} = 'V1 after; n=';
            elseif t == 3 % PSS before training
                if isempty(pssbf{a})
                    lbl{a,t,f} = '';
                    lblBit(t) = false;
                    continue
                end
                linStyl = '--';
                col = 'r';
                mrk = 'o';
                mrkSz = 5;
                if f == 1
                    tbl = pssbf{a};
                elseif f == 2
                    tbl = pssbf{a}(pssbf{a}.trainAx==1,:);
                elseif f == 3
                    tbl = pssbf{a}(pssbf{a}.orthAx==1,:);
                end
                lbl{a,t,f} = 'PSS before; n=';
            elseif t == 4 % PSS after training
                if isempty(pssaf{a})
                    lbl{a,t,f} = '';
                    lblBit(t) = false;
                    continue
                end
                linStyl = '-';
                col = 'r';
                mrk = '.';
                mrkSz = 20;
                if f == 1
                    tbl = pssaf{a};
                elseif f == 2
                    tbl = pssaf{a}(pssaf{a}.trainAx==1,:);
                elseif f == 3
                    tbl = pssaf{a}(pssaf{a}.orthAx==1,:);
                end
                lbl{a,t,f} = 'PSS after; n=';
            end

            tbl = tbl(tbl.goodUnit == 1,:);
            nU = height(tbl); 
            lbl{a,t,f} = [lbl{a,t,f} num2str(nU)];
            if nU == 0
                lblBit(t) = false;
                continue
            else
                lblBit(t) = true;
            end

            subplot(nAnimals+1,nMet+1,m+((nMet+1)*(a-1))); hold on

            if m<nMet+1

                metric = metrics{m};
                if strcmp(metric,'dsi')
                    dist{a,t,f}(:,m) = tbl.dsi;
                    xLbl = 'DSI';
                elseif strcmp(metric,'dcv')
                    dist{a,t,f}(:,m) = 1-tbl.dcv;
                    xLbl = '1-DCV';
                elseif strcmp(metric,'rPref')
                    dist{a,t,f}(:,m) = tbl.rPref;
                    xLbl = 'rPref';
                elseif strcmp(metric,'osi')
                    dist{a,t,f}(:,m) = tbl.osi;
                    xLbl = 'osi';
                end

                c(t) = cdfplot(dist{a,t,f}(:,m));
                c(t).LineStyle = linStyl;
                c(t).Color = col;
                c(t).LineWidth = 2;
                title(animals{a})
                xlabel(xLbl); 
                if strcmp(metric,'dsi') || strcmp(metric,'dcv') || strcmp(metric,'osi')
                    xlim([0 1])
                end
                ylabel('percentile')
            
                if t == 4 && m==1
                    legend(c(lblBit),lbl{a,lblBit,f},'Location','southeast')
                    clear c lblBit
                end

            elseif m == nMet+1
        
                x = mean(vertcat(tbl.tuningX{:}));
                y = mean(vertcat(tbl.tuningY{:}));
%             sem = std(vertcat(tbl.tuningY{:}))/sqrt(nU);
%             plot( x , y , [col linStyl] , 'LineWidth' , 2)
%             plot( repmat(x,2,1) , y+([-1;1]*sem) , col , 'LineWidth' , 2)
%             xlabel('ori (deg) relative to pref'); xticks([-180 -90 0 90 180])
%             ylabel('BCFR (Hz)')
%         
%             subplot(nAnimals+1,4,4+(4*(a-1))); hold on
                y = y/max(y);
                plot( x , y , [col linStyl] , 'LineWidth' , 2)
                xlabel('ori (deg) relative to pref'); xticks([-180 -90 0 90 180])
                ylabel('norm BCFR (Hz)')

            end
        
            clear tbl
        end
    end
    end
    
    sgtitle(ttl{f})
    set(gcf,'Position',[10 10 1010 1010])
    if svePlt == 1
        saveas(gcf,fullfile(figDir,['anaTrainDSI ' ttl{f}]),'fig')
    end

end


%% Stats

varNames = {'animal','units','metric','area/training','test1','h1','p1','test2','h2','p2'};
compCount = 0;
for a = 1:nAnimals+1
    for f = 1:nFig
        for m = 1:nMet
        
            compCount = compCount+1;
            tst1{compCount,1} = 'kstest2';
            tst2{compCount,1} = 'ranksum';
            if ~isempty(v1bf{a}) && ~isempty(dist{a,1,f})
                i = dist{a,1,f}(:,m);
            else
                i = [];
            end
            if ~isempty(v1af{a}) && ~isempty(dist{a,2,f})
                j = dist{a,2,f}(:,m);
            else
                j = [];
            end
            if ~isempty(i) && ~isempty(j)
                [h1(compCount),p1(compCount)] = kstest2(i,j);
                [p2(compCount),h2(compCount)] = ranksum(i,j);
            else
                h1(compCount) = 0;
                p1(compCount) = -1;
                h2(compCount) = 0;
                p2(compCount) = -1;
            end
            aID{compCount,1} = animals{a};
            fID{compCount,1} = ttl{f};
            mID{compCount,1} = metrics{m};
            tID{compCount,1} = 'v1bf VS v1af';
        
            compCount = compCount+1;
            tst1{compCount,1} = 'kstest2';
            tst2{compCount,1} = 'rankSum';
            if ~isempty(pssbf{a}) && ~isempty(dist{a,3,f})
                i = dist{a,3,f}(:,m);
            else
                i = [];
            end
            if ~isempty(pssaf{a}) && ~isempty(dist{a,4,f}) 
                j = dist{a,4,f}(:,m);
            else
                j = [];
            end
            if ~isempty(i) && ~isempty(j)
                [h1(compCount),p1(compCount)] = kstest2(i,j);
                [p2(compCount),h2(compCount)] = ranksum(i,j);
            else
                h1(compCount) = 0;
                p1(compCount) = -1;
                h2(compCount) = 0;
                p2(compCount) = -1;
            end
            aID{compCount,1} = animals{a};
            fID{compCount,1} = ttl{f};
            mID{compCount,1} = metrics{m};
            tID{compCount,1} = 'pssbf VS pssaf';


        end
    end
end
stats = table(aID,fID,mID,tID,tst1,h1',p1',tst2,h2',p2','VariableNames',varNames);
