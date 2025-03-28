clear all
close all

% dataFold = '/Volumes/Lab drive/Brandon/data/dataSets/DSdev';
% dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data/dataSets/DSdev';
% dataFold = 'F:\Brandon\data\dataSets\DSdev';
dataFold = 'Y:\Brandon\data\dataSets\DSdev';
load(fullfile(dataFold,"DSdev_dataSet.mat"))

areas = {'V1','PSS'};
ageGroups = {[29 32],[33 36],[37 300]};
% ageGroups = {[28 32],[29 33],[30 34],[31 35],[32 36],[33 37],[34 38],[35 39],[36 40],[37 41],[38 42],[39 43],[40 44],[41 300]};

nAG = length(ageGroups);
nAR = length(areas);
for ar = 1:nAR
for ag = 1:nAG
    ageLims = ageGroups{ag};
    areaIdx = strcmp(projectTbl.recSite,areas{ar});
    ageLimIdx = projectTbl.age>=ageLims(1) & projectTbl.age<=ageLims(2);
    dat{ar,ag} = vertcat(projectTbl(areaIdx & ageLimIdx,:).sumStats{:});

    dat{ar,ag} = dat{ar,ag}(dat{ar,ag}.goodUnit,:);

    [~,sortIdx] = sort(dat{ar,ag}.oriPref);
    dat{ar,ag} = dat{ar,ag}(sortIdx,:);
    nU(ar,ag) = height(dat{ar,ag});

    [cMean{ar,ag},cTrial{ar,ag},rMean{ar,ag},rTrial{ar,ag},~,~,score{ar,ag},...
        ~,D{ar,ag},Dshift{ar,ag},distF{ar,ag},distNull{ar,ag}] = anaPCA(dat{ar,ag});


    [~] = popDecode(rTrial{ar,ag},cTrial{ar,ag});

end
end

%plot

figure;
for ar = 1:nAR
    if strcmp(areas{ar},'V1')
        clr = 'b';
    elseif strcmp(areas{ar},'PSS')
        clr = 'r';
    end
    for ag = 1:nAG
    
        subplot(1,nAG)

    end
end

