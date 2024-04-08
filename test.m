clear all
close all

load('/Volumes/Lab drive/Brandon/data/dataSets/training/CoolV1Train_DPI_MUdataSet.mat')

figure
r = (1:8)/8;
colorsPSS = vertcat((1:8)/8,zeros(2,8));
colorsV1 = vertcat(zeros(2,8),(1:8)/8);
nAnimals = length(animals);
for a = 1:nAnimals
for h = 1:8
    V1 = vertcat(v1{1+(2*(h-1)),a} , v1{2+(2*(h-1)),a});
    V1 = V1(V1.goodUnit==1,:);
    if ~isempty(V1)
    subplot(nAnimals,2,1+(2*(a-1))); hold on
    cdf = cdfplot(V1.dpi);
    cdf.Color = colorsV1(:,h);
    cdf.LineWidth = 2;
    end

    PSS = vertcat(pss{1+(2*(h-1)),a} , pss{2+(2*(h-1)),a});
    PSS = PSS(PSS.goodUnit==1,:);
    if ~isempty(PSS)
    subplot(nAnimals,2,2+(2*(a-1))); hold on
    cdf = cdfplot(PSS.dpi);
    cdf.Color = colorsPSS(:,h);
    cdf.LineWidth = 2;
    end
    

%     subplot(2,1,1);hold on
%     [v] = violin(V1.dpi,0.2,[0 1 100],1,0);
%     plot((0.5*v.kdeX)+h,v.kdeY,'Color','b','LineWidth',2)
%     plot((-0.5*v.kdeX)+h,v.kdeY,'Color','b','LineWidth',2)
%     patch([0.5*v.kdeX -0.5*fliplr(v.kdeX)]+h,[v.kdeY fliplr(v.kdeY)],'b','EdgeColor','none','FaceAlpha',0.2)
%     plot([h-0.5 h+0.5],repmat(mean(V1.dpi),1,2),'k','LineWidth',2)
%     plot([h-0.5 h+0.5],repmat(median(V1.dpi),1,2),'Color',[0.8 0.8 0.8],'LineWidth',2)
% 
%     subplot(2,1,2);hold on
%     [v] = violin(PSS.dpi,0.2,[0 1 100],1,0);
%     plot((0.5*v.kdeX)+h,v.kdeY,'Color','r','LineWidth',2)
%     plot((-0.5*v.kdeX)+h,v.kdeY,'Color','r','LineWidth',2)
%     patch([0.5*v.kdeX -0.5*fliplr(v.kdeX)]+h,[v.kdeY fliplr(v.kdeY)],'r','EdgeColor','none','FaceAlpha',0.2)
%     plot([h-0.5 h+0.5],repmat(mean(PSS.dpi),1,2),'k','LineWidth',2)
%     plot([h-0.5 h+0.5],repmat(median(PSS.dpi),1,2),'Color',[0.8 0.8 0.8],'LineWidth',2)

end
end