%population vector

clear all
close all

anaMode = 'MU';
trainProj = 'Train_V1Cool';
% dataFold = 'F:\Brandon\data';
dataFold = '/Volumes/Lab drive/Brandon/data';
dsFold = fullfile(dataFold,'dataSets','training',trainProj,anaMode);
load(fullfile(dsFold,'ranksum & rPref above 2',[trainProj '_' anaMode 'dataSet.mat']))

stim = 180;
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
    
        angDir(u) = tbl.meanVec{u}.angDir;
        ldr(u) = tbl.meanVec{u}.ldr;
        oriInd = find(strcmp(tbl.paramKey{u},'ori'));
        tune(u,:) = mean(tbl.response{u},'omitnan');
        ori(u,:) = tbl.condition{u}(oriInd,:);
        
    
    %     polarplot(repmat(deg2rad(angDir),1,2),[0 ldr],'k');hold on
    
    end
    
    w = tune(ori==stim);
    pV = sum(w.*exp(sqrt(-1)*mod(deg2rad(angDir(:)),2*pi)));
    
    
    PV(i) = rad2deg(mod( angle(pV) ,2*pi));


end
