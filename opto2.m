clear all
close all

anaMode = 'MU';
dataFold = 'F:\Brandon\data';
animal = 'febs6'; unit = '000'; expt = '000'; probe = 1;
baseName = [animal '_u' unit '_' expt];

physDir = fullfile(dataFold,'Ephys',animal,baseName);
load(fullfile(physDir,[baseName '_trialInfo.mat']))
load(fullfile(physDir,[baseName '.analyzer']),'-mat','Analyzer')
load(fullfile(physDir,[baseName '_id.mat']))

predelay = getparam('predelay',Analyzer);
stimTime = getparam('stim_time',Analyzer);
postdelay = getparam('postdelay',Analyzer);
stimStartID = 3;
tLims = [1,stimTime];

switch anaMode
    case 'SU'
SUTrialData(fullfile(dataFold,'Ephys'),animal,unit,expt,probe,'id',stimStartID,tLims(1),tLims(2),0)
load(fullfile(physDir,[baseName '_p' num2str(probe) '_SUTrial.mat']),'SU','SUinfo')
uDat = SU; uInfo = SUinfo;
clear SU SUinfo
    case 'MU'
MUThreshTrialData(fullfile(dataFold,'Ephys'),animal,unit,expt,probe,'id',stimStartID,tLims(1),tLims(2),0)
load(fullfile(physDir,[baseName '_p' num2str(probe) '_MUThreshTrial.mat']),'MUThresh','MUThreshInfo')
uDat = MUThresh; uInfo = MUThreshInfo;
clear MUThresh MUThreshInfo
end
