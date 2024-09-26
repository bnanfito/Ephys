%anaFI

clear all
close all

anaMode = 'MU';
trainProj = 'Train_V1Cool';
% dataFold = 'F:\Brandon\data';
dataFold = '/Volumes/Lab drive/Brandon/data';
dsFold = fullfile(dataFold,'dataSets','training',trainProj,anaMode);
load(fullfile(dsFold,[trainProj '_' anaMode 'dataSet.mat']))

figure;hold on
for i = 1:4

    if i==1
        I = 'v1bf';
        tbl = data.v1bf;
        clr = 'b';
        linStyl = '--';
    elseif i==2
        I = 'v1af';
        tbl = data.v1af;
        clr = 'b';
        linStyl = '-';
    elseif i==3
        I = 'pssbf';
        tbl = data.pssbf;
        clr = 'r';
        linStyl = '--';
    elseif i==4
        I = 'pssaf';
        tbl = data.pssaf;
        clr = 'r';
        linStyl = '-';
    end

    nU = height(tbl);
    for u = 1:nU

        r = tbl.response{u};
        c = tbl.condition{u};

        FI{i}(u) = fisherInfo(r,c);

    end

    FI{i} = FI{i}(tbl.goodUnit);
    cdf = cdfplot(FI{i});
    cdf.LineStyle = linStyl;
    cdf.Color = clr;
    cdf.LineWidth = 2;


end