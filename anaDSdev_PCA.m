clear all
close all

% dataFold = '/Volumes/Lab drive/Brandon/data/dataSets/DSdev';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/DSdev';
% dataFold = 'F:\Brandon\data\dataSets\DSdev';
dataFold = 'Y:\Brandon\data\dataSets\DSdev';
load(fullfile(dataFold,"DSdev_dataSet.mat"))

area = 'V1';
ageGroups = {[29 32],[33 36],[37 300]};
% ageGroups = {[28 32],[29 33],[30 34],[31 35],[32 36],[33 37],[34 38],[35 39],[36 40],[37 41],[38 42],[39 43],[40 44],[41 300]};

nAG = length(ageGroups);
for ag = 1
    ageLims = ageGroups{ag};
    areaIdx = strcmp(projectTbl.recSite,area);
    ageLimIdx = projectTbl.age>=ageLims(1) & projectTbl.age<=ageLims(2);
    dat{ag} = vertcat(projectTbl(areaIdx & ageLimIdx,:).sumStats{:});

    dat{ag} = dat{ag}(dat{ag}.goodUnit,:);

    [~,sortIdx] = sort(dat{ag}.oriPref);
    dat{ag} = dat{ag}(sortIdx,:);
    nU(ag) = height(dat{ag});

    [cMean{ag},cTrial{ag},rMean{ag},rTrial{ag},~,~,score{ag},~,D{1,ag},Dshift{1,ag},distF{1,ag},distNull{1,ag}] = anaPCA(dat{ag});
    rTrial{ag} = rTrial{ag}./max(rTrial{ag});

    [~] = popDecode(rTrial{ag},cTrial{ag});

end


