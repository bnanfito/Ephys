function anaExpt(animal,unit,expt,probe,anaMode,dataFold,plotBit)

exptName = [animal '_u' unit '_' expt];
physDir = fullfile(dataFold,'Ephys',animal,exptName);
load(fullfile(physDir,[exptName '_id.mat']))
if strcmp(anaMode,'MU') 
    fileName = fullfile(physDir,[exptName '_p' num2str(probe) '_MUThreshTrial.mat']);
    if ~isfile(fileName)
        MUThreshTrialData(fullfile(dataFold,'Ephys'),animal,unit,expt,probe,'id',3,1,1)
    end
    load(fileName,'MUThresh','MUThreshInfo')
    data = MUThresh;
    info = MUThreshInfo;
    clear MUThresh MUThreshInfo
elseif strcmp(anaMode,'SU')
    fileName = fullfile(physDir,[exptName '_p' num2str(probe) '_SUTrial.mat']);
    if ~isfile(fileName)
        SUTrialData(fullfile(dataFold,'Ephys'),animal,unit,expt,probe,'id',3,1,1)
    end
    load(fileName,'SU','SUinfo')
    data = SU;
    info = SUinfo;
    clear SU SUinfo
end

nUnits = length(data);
for u = nUnits

end




end