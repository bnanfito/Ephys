%plotAnimals
clear all
close all

anaMode = 'SU';

%% exclude bi ff

badAnimals = {'febe7','febg6','febg7','febg8','febg9','febk6','febl6','febl8','febi6','febg3','febe2','febl4','febg0','febf7','febg1','febg2','febf6','febe0'};
badExpts = {'febh5_u001_012'};
biFFlist = {'febk7_u000_002','febk8_u000_000','febl0_u000_001','febh5_u000_000','febj8_u000_000','febj8_u001_000','febj8_u002_000','febj8_u003_000'};

needSort = {'febj8_u003_009','febj8_u003_010','febj8_u003_022','febl0_u000_012','febl0_u000_023','febl0_u000_027','febl0_u001_010','febl0_u001_019','febj7_u000_011'};

%% query data

[coolExpts] = queryCoolV1([0 500]);
for f = 1:height(coolExpts)
    exptNameList{f,1} = [coolExpts.experimentId{f} '_u' coolExpts.unitNr{f} '_' coolExpts.experimentNr{f}];
end
coolExpts.exptName = exptNameList;
excludeIdx = ismember(coolExpts.experimentId,badAnimals) |...
                ismember(coolExpts.exptName,badExpts) |...
                ismember(coolExpts.exptName,biFFlist) |...
                ismember(coolExpts.exptName,needSort);
pssIdx = strcmp(coolExpts.recSite,'PSS');
ori16Idx = strcmp(coolExpts.module,'PG') &...
       (strcmp(coolExpts.looperNameCond1,'ori16') |...
        strcmp(coolExpts.looperNameCond2,'ori16') |...
        strcmp(coolExpts.looperNameCond3,'ori16'));
coolExpts = coolExpts(~excludeIdx & pssIdx & ori16Idx,:);
animals = unique(coolExpts.experimentId);
for a = 1:length(animals)
    ages(a,1) = unique( coolExpts.age( strcmp(coolExpts.experimentId,animals(a)) ) );
end

clear exptNameList

[cntrlExpts] = queryCntrl;
for f = 1:height(cntrlExpts)
    exptNameList{f,1} = [cntrlExpts.experimentId{f} '_u' cntrlExpts.unitNr{f} '_' cntrlExpts.experimentNr{f}];
end
cntrlExpts.exptName = exptNameList;
includeIdx = ismember(cntrlExpts.experimentId,animals);
excludeIdx = ismember(cntrlExpts.exptName,badExpts) |...
                ismember(cntrlExpts.exptName,biFFlist) |...
                ismember(cntrlExpts.exptName,needSort);
pssIdx = strcmp(cntrlExpts.recSite,'PSS');
ori16Idx = strcmp(cntrlExpts.module,'PG') &...
       (strcmp(cntrlExpts.looperNameCond1,'ori16') |...
        strcmp(cntrlExpts.looperNameCond2,'ori16') |...
        strcmp(cntrlExpts.looperNameCond3,'ori16'));
cntrlExpts = cntrlExpts(~excludeIdx & includeIdx & pssIdx & ori16Idx,:);

clear exptNameList

%% coolExpts

for f = 1:height(coolExpts)

    disp(['generating sumStats for: ' coolExpts.experimentId{f} '_u' coolExpts.unitNr{f} '_' coolExpts.experimentNr{f} '_p' num2str(coolExpts.probeId(f))])
    [sumStats] = anaOri(coolExpts.experimentId{f},...
                        coolExpts.unitNr{f},...
                        coolExpts.experimentNr{f},...
                        coolExpts.probeId(f),...
                        anaMode,'mono c hemi',0,0);
    sumStats = sumStats(sumStats.goodUnit,:);
    RP{f,1} = sumStats.rPref;
    LDir{f,1} = sumStats.ldr;

    clear sumStats

end

coolExpts.RP = RP;
coolExpts.LDir = LDir;

clear RP LDir 

%% cntrlExpts

for f = 1:height(cntrlExpts)

    disp(['generating sumStats for: ' cntrlExpts.experimentId{f} '_u' cntrlExpts.unitNr{f} '_' cntrlExpts.experimentNr{f} '_p' num2str(cntrlExpts.probeId(f))])
    [sumStats] = anaOri(cntrlExpts.experimentId{f},...
                        cntrlExpts.unitNr{f},...
                        cntrlExpts.experimentNr{f},...
                        cntrlExpts.probeId(f),...
                        anaMode,'mono c hemi',0,0);
    sumStats = sumStats(sumStats.goodUnit,:);
    RP{f,1} = sumStats.rPref;
    LDir{f,1} = sumStats.ldr;

    clear sumStats

end

cntrlExpts.RP = RP;
cntrlExpts.LDir = LDir;

clear RP LDir

%% animal averages


for a = 1:length(animals)

    curAni = animals{a};
    
    aniIdxCool = strcmp(coolExpts.experimentId,curAni);
    aniCoolExpts = coolExpts(aniIdxCool,:);
    aniIdxCntrl = strcmp(cntrlExpts.experimentId,curAni);
    aniCntrlExpts = cntrlExpts(aniIdxCntrl,:);

    pens = unique(aniCoolExpts.penNr)';
    coolRP = [];
    coolLdir = [];
    for p = pens
        curPenExpts = aniCoolExpts(aniCoolExpts.penNr==p,:);
        for f = 1:height(curPenExpts)
            if f == 1
                fIdx = f;
            elseif length(curPenExpts.RP{f})>length(curPenExpts.RP{1})
                fIdx = f;
            end
        end
        coolRP = [coolRP;curPenExpts.RP{fIdx}];
        coolLdir = [coolLdir;curPenExpts.LDir{fIdx}];
        aniExptsCool{a,p} = curPenExpts.exptName{fIdx};
        clear curPenExpts
    end
    aniRPcool(a) = mean(coolRP,'omitnan');
    aniLdirCool(a) = mean(coolLdir,'omitnan');
    nCool(a) = length(coolRP);
    clear coolRP coolLdir

    pens = unique(aniCntrlExpts.penNr)';
    cntrlRP = [];
    cntrlLdir = [];
    for p = pens
        curPenExpts = aniCntrlExpts(aniCntrlExpts.penNr==p,:);
        for f = 1:height(curPenExpts)
            if f == 1
                fIdx = f;
            elseif length(curPenExpts.RP{f})>length(curPenExpts.RP{1})
                fIdx = f;
            end
        end
        cntrlRP = [cntrlRP;curPenExpts.RP{fIdx}];
        cntrlLdir = [cntrlLdir;curPenExpts.LDir{fIdx}];
        aniExptsCntrl{a,p} = curPenExpts.exptName{fIdx};
        clear curPenExpts
    end
    aniRPcntrl(a) = mean(cntrlRP,'omitnan');
    aniLdirCntrl(a) = mean(cntrlLdir,'omitnan');
    nCntrl(a) = length(cntrlRP);
    clear cntrlRP cntrlLdir
end

%% Plot

if strcmp(anaMode,'MU')

    figure; hold on
    for a = 1:length(animals)
        plot(ages(a),aniRPcntrl(a),'k.','MarkerSize',nCntrl(a))
        plot(ages(a),aniRPcool(a),'c.','MarkerSize',nCool(a))
    end
    xlabel('age')
    ylabel('Rpref')
    
    
    figure; hold on
    for a = 1:length(animals)
        si = (aniRPcool(a) - aniRPcntrl(a))/(aniRPcool(a) + aniRPcntrl(a));
        plot(ages(a),si,'g.','MarkerSize',10)
    end
    xlabel('age')
    ylabel('SI')
    
    
    figure; hold on
    for a = 1:length(animals)
        plot(ages(a),aniLdirCntrl(a),'k.','MarkerSize',nCntrl(a))
        plot(ages(a),aniLdirCool(a),'c.','MarkerSize',nCool(a))
    end
    xlabel('age')
    ylabel('Ldir')

elseif strcmp(anaMode,'SU')

    mrkSize = 20;

    figure; hold on
    plot(ages,aniRPcntrl,'k.','MarkerSize',mrkSize)
    plot(ages,aniRPcool,'c.','MarkerSize',mrkSize)
    xlabel('age')
    ylabel('Rpref')
    
    
    figure; hold on
    si = (aniRPcool - aniRPcntrl)./(aniRPcool + aniRPcntrl);
    plot(ages,si,'k.','MarkerSize',mrkSize)
    xlabel('age')
    ylabel('SI')
    
    
    figure; hold on
    plot(ages,aniLdirCntrl,'k.','MarkerSize',mrkSize)
    plot(ages,aniLdirCool,'c.','MarkerSize',mrkSize)
    xlabel('age')
    ylabel('Ldir')

end


