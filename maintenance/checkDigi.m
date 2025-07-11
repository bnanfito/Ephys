clear all
close all
exptFold = 'Z:\EphysNew\data\feav2';
cd(exptFold)
expts = dir;
expts = expts(vertcat(expts.isdir)==1,:);
expts = expts(3:end,:);

countGood = 0;
countBad = 0;
for e = 1:size(expts,1)

    baseName = expts(e).name;
    disp(baseName)
    cd(fullfile(exptFold,baseName))

    physname = fullfile(exptFold,baseName,[baseName '_digitalin.dat']);
    DigiFile = fopen(physname);
    fileinfo = dir(physname);
    digiData = fread(DigiFile, (fileinfo.bytes)/2, 'uint16');
    nDigi(e,1) = length(unique(digiData));
    if nDigi(e,1) > 2
        countGood = countGood+1;
        goodExpts{countGood,1} = baseName;
    elseif nDigi(e,1) < 2
        countBad = countBad+1;
        badExpts{countBad,1} = baseName;
    end
    fclose(DigiFile);

end