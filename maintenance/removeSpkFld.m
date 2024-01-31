
animal = 'FEAM7';
unit = '000';
exp = '001';
probe = 1;
nJobs = 100;
fileBase = [animal '_u' unit '_' exp];
dataFold = '/Users/brandonnanfito/OneDrive - Johns Hopkins/Documents/data';
field = {'PC1','PC2','PC3'};

for job = 0:nJobs-1
    disp(num2str(job))
    
    load(fullfile(dataFold,animal,fileBase,'SpikeFiles',[fileBase '_j' num2str(job) '_p' num2str(probe) '_spkinfo.mat']))
    if isfield(spk,field)
        spk = rmfield(spk,field);
        save([fileBase '_j' num2str(job) '_p' num2str(probe) '_spkinfo'],'spk')
    else
        disp(['job' num2str(job) 'doesnt have this field'])
    end
    
end
