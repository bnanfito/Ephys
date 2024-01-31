
animal = 'febe9';
unit = '001';
expt = '008';
baseName = [animal '_u' unit '_' expt];
if ispc
    dataFold = 'C:\Users\brand\OneDrive - Johns Hopkins\Documents\data\Ephys';
elseif ismac
    dataFold = '/Users/brandonnanfito/Library/CloudStorage/OneDrive-JohnsHopkins/Documents/data/Ephys';
end
cd(fullfile(dataFold,animal,baseName))
load(fullfile(dataFold,animal,baseName,[baseName '_id.mat']))

for p =1:length(id.probes)
    if ~isempty(id.probes(p).type)
        probewiring=feval(['probeConfig_' id.probes(p).type]);
        id.probes(p).channels = probewiring(:,1);
        id.probes(p).x = probewiring(:,2);
        id.probes(p).y = probewiring(:,3);
        id.probes(p).z = probewiring(:,4);
        id.probes(p).shaft = probewiring(:,5);
        id.probes(p).nChannels = length(probewiring(:,1));
    end
end
save([baseName '_id.mat'],'id')