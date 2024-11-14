clear all
close all

animal = 'febi0';
unit = {'000','000'};
expt = {'007','008'};
area = 'V1';
anaMode = 'MU';
dataFold = 'F:\Brandon\data';

figure;hold on
for e = 1:length(expt)

    if strcmp(expt{e},'007')
        clr = 'k';
        linStyl = '-';
        legLbl{e} = 'R eye';
    elseif strcmp(expt{e},'008')
        clr = 'k';
        linStyl = '--';
        legLbl{e} = 'L eye';
    end
    
    [sumStats{e}] = anaSF(animal,unit{e},expt{e},area,anaMode,dataFold,0,0);
    sumStats{e} = sumStats{e}(sumStats{e}.goodUnit,:);

    nU = height(sumStats{e});
    for u = 1:nU

        oriInd = strcmp(sumStats{e}.paramKey{u},'ori');
        sfInd = strcmp(sumStats{e}.paramKey{u},'s_freq');

        opIdx = squeeze(mean(sumStats{e}.condition{u}(oriInd,:,:),2)) == sumStats{e}.oriPref(u);
        y(u,:) = mean( sumStats{e}.response{u}(:,:,opIdx) ,1,'omitnan');
        y(u,:) = y(u,:)./max(y(u,:));
        x = sumStats{e}.condition{u}(sfInd,:,opIdx);

    end

    p(e) = plot(x,mean(y,'omitnan'),[linStyl clr],'LineWidth',2);
    sem = std(y,'omitnan')/sqrt(size(y,1));
    plot(repmat(x,2,1),mean(y)+([1;-1].*sem),clr,'LineWidth',2)

end

legend(p,legLbl)
