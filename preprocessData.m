
function preprocessData(animal,unit,expt,probe,anaMode,copyToZ,dataFold)

if strcmp(anaMode,'SU')
    MUflag = 0;
elseif strcmp(anaMode,'MU')
    MUflag = 1;
end

physDir = fullfile(dataFold,'Ephys');
extractTrials(physDir,physDir,physDir,animal,unit,expt)
name = 'BRN';
scaleFactor = 4;
nJobs = 100;
legacyFlag = 0;

baseName = [animal '_u' unit '_' expt];

load(fullfile(physDir,animal,baseName,[baseName '_id.mat']),'id');
nChannels=sum([id.probes.nChannels]);

hName=fullfile(physDir,animal,baseName,[baseName '_info.rhd']);
h=read_Intan_Header(hName);
sf = h.sample_rate;
fileinfo = dir(fullfile(physDir,animal,baseName,[baseName '_amplifier.dat']));
samples = fileinfo.bytes/(2*nChannels); % Number of samples in amplifier data file

for p = probe

    disp(['processing ' baseName ' probe ' num2str(p)])

    % Compute thresholds automatically if MUA; otherwise for SUA threshold
    % file should be generated with thresholdGui and should be found in
    % physDir > animal > baseName
    if MUflag == 1
        computeMUThreshold(physDir,animal,unit,expt,p,(samples/sf)*.25,scaleFactor,name,copyToZ)
    end

    % spike detection
    parfor j = 0:nJobs-1
        extractSpikes(physDir,animal,unit,expt,p,name,copyToZ,MUflag,legacyFlag,nJobs,j)
    end

    % feature extraction
    parfor j = 0:nJobs-1
        extractSpikeProps(physDir,animal,unit,expt,p,name,copyToZ,MUflag,j)
    end

    if MUflag == 1 % for MUA generate MUspkMerge.mat file automatically, otherwise for SUA generate spkSort.mat file using sortGui (manual)
        mergeMUspkInfo(physDir,animal,unit,expt,p,0,nJobs-1)
    else % extra feature extraction for SUA spk sorting: compute PCs based on waveforms
        wavePCA(animal,unit,expt,p,nJobs,0.5,physDir)
        parfor j = 0:nJobs-1
            spkWvPC(animal,unit,expt,p,j,physDir)
        end
    end



end

fclose all;
end

