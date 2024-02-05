

function [dataRaw,dataFilt] = readRawData(animal,unit,expt,start,stop,band)

if ispc
%     dataFold = 'C:\Users\brand\Documents\data\Ephys';
    dataFold = 'D:\OneDrive - Johns Hopkins\Documents\data\Ephys';
elseif ismac
    dataFold = '/Volumes/Lab drive/Brandon/Ephys';
end
baseName = [animal '_u' unit '_' expt];
cd(fullfile(dataFold,animal,baseName))

load(fullfile(dataFold,animal,baseName,[baseName '_id.mat']),'id')

nProbes = length(id.probes);
nCh = nan(nProbes,1);
for p = 1:nProbes
    nCh(p) = id.probes(p).nChannels;
end
sf= id.sampleFreq;


ampFile = [baseName '_amplifier.dat'];

startSample = start-1;
nSamps = (stop-start)+1;

fid = fopen(ampFile,'r');
fseek(fid,2*sum(nCh)*startSample,'bof');
dataRaw = fread(fid, [sum(nCh) nSamps], 'int16');

hp=band(1);
lp=band(2);
[b1,a1]=butter(3,[hp/sf,lp/sf]*2,'bandpass');
dataFilt=filter(b1,a1,dataRaw');dataFilt = dataFilt';


end