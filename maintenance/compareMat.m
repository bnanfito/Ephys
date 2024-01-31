%compare mat files
function out = compareMat(animal,unit,expt,probe)

% animal = 'FEAN6';
% unit = '000';
% expt = '000';
% probe = 1;
exptName = [animal '_u' unit '_' expt];
pathA = fullfile('F:\Brandon\data\Ephys',animal,exptName);
pathB = fullfile('D:\data\Ephys',animal,exptName);

a = load(fullfile(pathA,[exptName '_p' num2str(probe) '_MUspkMerge.mat']));
b = load(fullfile(pathB,[exptName '_p' num2str(probe) '_MUspkMerge.mat']));

out = unique(a.MUspkMerge.spktimes == b.MUspkMerge.spktimes);

end