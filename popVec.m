%Population Vector

clear all
close all

anaMode = 'SU';
load(['/Volumes/Lab drive/Brandon/data/dataSets/training/Train_V1Cool/' anaMode '/Train_V1Cool_' anaMode 'dataSet.mat'])


figure;hold on
for i = 1

    if i==1
        tbl = data.v1bf;
        clr = 'b';
        linStyl = '--';
    elseif i==2
        tbl = data.v1af;
        clr = 'b';
        linStyl = '-';
    elseif i==3
        tbl = data.pssbf;
        clr = 'r';
        linStyl = '--';
    elseif i==4
        tbl = data.pssaf;
        clr = 'r';
        linStyl = '-';
    end







end