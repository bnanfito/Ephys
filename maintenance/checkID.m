%Check ID
clear all
close all

if ispc
    dataFold = 'D:\OneDrive - Johns Hopkins\Documents\data';
%     dataFold = 'F:\Brandon\data\sf4rs01';
elseif ismac
    dataFold = '/Users/brandonnanfito/Library/CloudStorage/OneDrive-JohnsHopkins/Documents/data';
end
physDir = fullfile(dataFold,'Ephys');

% animals = {'FEAO4','FEAN6','FEAS6','FEAT1','FEAQ5'};
% animals = {'febj4'};
animals = {'febj3'};
nAnimals = length(animals);

for a = 1:nAnimals
    animal = animals{a};
%     exptList = dir(fullfile(physDir,animal));

    if strcmp(animal,'FEAO4')
    
        blocks{1,a} = 'FEAO4_u001_000';
        blocks{2,a} = 'FEAO4_u001_001';
        blocks{3,a} = 'FEAO4_u001_002';
        blocks{4,a} = 'FEAO4_u001_003';
        blocks{5,a} = 'FEAO4_u001_004';
        blocks{6,a} = 'FEAO4_u001_005';
        blocks{7,a} = 'FEAO4_u001_006';
        blocks{8,a} = 'FEAO4_u001_007';
        blocks{9,a} = 'FEAO4_u001_008';
        blocks{10,a} = 'FEAO4_u001_009';
        blocks{11,a} = 'FEAO4_u001_010';
        blocks{12,a} = 'FEAO4_u001_011';
        blocks{13,a} = 'FEAO4_u001_012';
        blocks{14,a} = 'FEAO4_u001_013';
        blocks{15,a} = 'FEAO4_u001_015';
        blocks{16,a} = 'FEAO4_u001_016';
    
        flagBlocks = [];
    
    elseif strcmp(animal,'FEAN6')
    
        blocks{1,a} = 'FEAN6_u001_000';
        blocks{2,a} = 'FEAN6_u001_001';
        blocks{3,a} = 'FEAN6_u001_002';
        blocks{4,a} = 'FEAN6_u001_003';
        blocks{5,a} = 'FEAN6_u001_004';
        blocks{6,a} = 'FEAN6_u001_005';
        blocks{7,a} = 'FEAN6_u001_006';
        blocks{8,a} = 'FEAN6_u001_007';
        blocks{9,a} = 'FEAN6_u001_008';
        blocks{10,a} = 'FEAN6_u001_009';
        blocks{11,a} = 'FEAN6_u001_010';
        blocks{12,a} = 'FEAN6_u001_011';
        blocks{13,a} = 'FEAN6_u001_012';
        blocks{14,a} = 'FEAN6_u001_013';
        blocks{15,a} = 'FEAN6_u001_014';
        blocks{16,a} = 'FEAN6_u001_015';
    
        flagBlocks = [];
    
    elseif strcmp(animal,'FEAS6')
    
        blocks{1,a} = 'FEAS6_u001_000';
        blocks{2,a} = 'FEAS6_u001_001';
        blocks{3,a} = 'FEAS6_u001_002';
        blocks{4,a} = 'FEAS6_u001_003';
        blocks{5,a} = 'FEAS6_u001_004';
        blocks{6,a} = 'FEAS6_u001_005';
        blocks{7,a} = 'FEAS6_u001_006';
        blocks{8,a} = 'FEAS6_u001_007';
        blocks{9,a} = 'FEAS6_u001_008';
        blocks{10,a} = 'FEAS6_u001_009'; %digi inputs missing
        blocks{11,a} = 'FEAS6_u001_010'; %digi inputs missing
        blocks{12,a} = 'FEAS6_u001_011'; %digi inputs missing
        blocks{13,a} = 'FEAS6_u001_012'; %digi inputs missing
        blocks{14,a} = 'FEAS6_u001_013'; %digi inputs missing
        blocks{15,a} = 'FEAS6_u001_014';
        blocks{16,a} = 'FEAS6_u001_015';
    
        flagBlocks = 10:14;
    
    elseif strcmp(animal,'FEAT1')
    
        blocks{1,a} = 'FEAT1_u001_000';
        blocks{2,a} = 'FEAT1_u001_001';
        blocks{3,a} = 'FEAT1_u001_002';
        blocks{4,a} = 'FEAT1_u001_003';
        blocks{5,a} = 'FEAT1_u001_004';
        blocks{6,a} = 'FEAT1_u001_005';
        blocks{7,a} = 'FEAT1_u001_006';
        blocks{8,a} = 'FEAT1_u001_007';
        blocks{9,a} = 'FEAT1_u001_008';
        blocks{10,a} = 'FEAT1_u001_009';
        blocks{11,a} = 'FEAT1_u001_010';
        blocks{12,a} = 'FEAT1_u001_011';
        blocks{13,a} = 'FEAT1_u001_012';
        blocks{14,a} = 'FEAT1_u001_013';
        blocks{15,a} = 'FEAT1_u001_014';
        blocks{16,a} = 'FEAT1_u001_015';
    
        flagBlocks = [];
    
    elseif strcmp(animal,'FEAQ5')
    
        blocks{1,a} = 'FEAQ5_u001_000';
        blocks{2,a} = 'FEAQ5_u001_001';
        blocks{3,a} = 'FEAQ5_u001_002';
        blocks{4,a} = 'FEAQ5_u001_003';
        blocks{5,a} = 'FEAQ5_u001_004';
        blocks{6,a} = 'FEAQ5_u001_005';
        blocks{7,a} = 'FEAQ5_u001_006';
        blocks{8,a} = 'FEAQ5_u001_007';
        blocks{9,a} = 'FEAQ5_u001_008';
        blocks{10,a} = 'FEAQ5_u001_009';
        blocks{11,a} = 'FEAQ5_u001_010';
        blocks{12,a} = 'FEAQ5_u001_011';
        blocks{13,a} = 'FEAQ5_u001_012';
        blocks{14,a} = 'FEAQ5_u001_013';
        blocks{15,a} = 'FEAQ5_u001_014';
        blocks{16,a} = 'FEAQ5_u001_015';
    
        flagBlocks = [];
    
    elseif strcmp(animal,'febh2')

        blocks{1,a} = 'febh2_u000_005'; %electrical noise issues (from pump)
        blocks{2,a} = 'febh2_u000_007'; %electrical noise issues (from pump)
        blocks{3,a} = 'febh2_u000_009';
        blocks{4,a} = 'febh2_u000_011';
        blocks{5,a} = 'febh2_u000_013';
        blocks{6,a} = 'febh2_u000_016';
        blocks{7,a} = 'febh2_u000_018'; %digi inputs missing
        blocks{8,a} = 'febh2_u000_020'; %digi inputs missing
        blocks{9,a} = 'febh2_u000_022'; %digi inputs missing
        blocks{10,a} = 'febh2_u000_024'; %digi inputs missing
        blocks{11,a} = 'febh2_u000_026'; %digi inputs missing
        blocks{12,a} = 'febh2_u000_028';
        blocks{13,a} = 'febh2_u000_030';
        blocks{14,a} = 'febh2_u000_032';
        blocks{15,a} = 'febh2_u000_034';
        blocks{16,a} = 'febh2_u000_036';

        flagBlocks = [];

    elseif strcmp(animal,'febh3')

        blocks{1,a} = 'febh3_u000_006';
        blocks{2,a} = 'febh3_u000_008';
        blocks{3,a} = 'febh3_u000_010';
        blocks{4,a} = 'febh3_u000_014'; %something weird with the trialInfo file
        blocks{5,a} = 'febh3_u000_017';
        blocks{6,a} = 'febh3_u000_019';
        blocks{7,a} = 'febh3_u000_021';
        blocks{8,a} = 'febh3_u000_023';
        blocks{9,a} = 'febh3_u000_025';
        blocks{10,a} = 'febh3_u000_027';
        blocks{11,a} = 'febh3_u000_029';
        blocks{12,a} = 'febh3_u000_031';
        blocks{13,a} = 'febh3_u000_033';
        blocks{14,a} = 'febh3_u000_035';
        blocks{15,a} = 'febh3_u000_037';
        blocks{16,a} = 'febh3_u000_039';

        flagBlocks = 4;

    elseif strcmp(animal,'febj3')

        blocks{1,a} = 'febj3_u000_005';
        blocks{2,a} = 'febj3_u000_007';
        blocks{3,a} = 'febj3_u000_009';
        blocks{4,a} = 'febj3_u000_011';
        blocks{5,a} = 'febj3_u000_013';
        blocks{6,a} = 'febj3_u000_015';
        blocks{7,a} = 'febj3_u000_017';
        blocks{8,a} = 'febj3_u000_019';
        blocks{9,a} = 'febj3_u000_021';
        blocks{10,a} = 'febj3_u000_023';
        blocks{11,a} = 'febj3_u000_025';
        blocks{12,a} = 'febj3_u000_028'; %fewer reps, block split with febj3_u000_027
        blocks{13,a} = 'febj3_u000_030';
        blocks{14,a} = 'febj3_u000_032';
        blocks{15,a} = 'febj3_u000_034';
        blocks{16,a} = 'febj3_u000_036';

        flagBlocks = [];

    elseif strcmp(animal,'febj4')

        blocks{1,a} = 'febj4_u000_008';
        blocks{2,a} = 'febj4_u000_009';
        blocks{3,a} = 'febj4_u000_010';
        blocks{4,a} = 'febj4_u000_011';
        blocks{5,a} = 'febj4_u000_012';
        blocks{6,a} = 'febj4_u000_013';
        blocks{7,a} = 'febj4_u000_014';
        blocks{8,a} = 'febj4_u000_015';
        blocks{9,a} = 'febj4_u000_016';
        blocks{10,a} = 'febj4_u000_017';
        blocks{11,a} = 'febj4_u000_018';
        blocks{12,a} = 'febj4_u000_019';
        blocks{13,a} = 'febj4_u000_020';
        blocks{14,a} = 'febj4_u000_021';
        blocks{15,a} = 'febj4_u000_022';
        blocks{16,a} = 'febj4_u000_023';

        flagBlocks = [];

    end

    for b = 1:size(blocks,1)

        exptName = blocks{b,a};
        if (~(contains(exptName,'fe') || contains(exptName,'FE')) || length(exptName)~=14) || ismember(b,flagBlocks)
            continue
        end
        disp(exptName)

        load(fullfile(physDir,animal,exptName,[exptName '_id.mat']))
        exptNames{b,a} = exptName;
        hasThreshSet(b,a) = isfield(id.MUthreshold,'settings');
        if isfield(id.MUthreshold,'settings')
        for p = 1:length(id.probes)
            threshLvl(b,a,p) = id.MUthreshold.settings{p}.scaleFactor;
        end
        end

    end

end



