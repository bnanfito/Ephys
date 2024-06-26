clear all; close all

area = 'PSS';
anaMode = 'MU';

if ispc
%     dataFold = 'D:\data'; 
%     dataFold = 'C:\Users\brand\Documents\data';
%     dataFold = 'F:\Brandon\data';
    dataFold = 'F:\Brandon\VSS2024\data';
elseif ismac
%     dataFold = '/Volumes/Lab drive/Brandon/data';
    dataFold = '/Volumes/Lab drive/Brandon/VSS2024/data';
%     dataFold = '/Users/brandonnanfito/Documents/NielsenLab/data';
end
physDir = fullfile(dataFold,'Ephys');
condNames = {'precool','before','during','after'};
condColor = {'k','k','c','m'};
conds = [1,3];

data(1).age = 'P30-32';
data(1).expList{1,1} =    {'febe7_u000_001'};
data(1).expList{2,1} =    {'febe7_u000_001'};
data(1).expList{3,1} =    {'febe7_u000_003'};
data(1).expList{4,1} =    {'febe7_u000_005'};

data(2).age = 'P33-35';
data(2).expList{1,1} =    {'febg7_u000_003';                    'febg8_u001_003'};
data(2).expList{2,1} =    {'febg7_u000_003';'febg7_u001_004';   'febg8_u001_003';'febg8_u002_002'};
data(2).expList{3,1} =    {'febg7_u000_008';'febg7_u001_005';   'febg8_u001_004';'febg8_u002_003'};
data(2).expList{4,1} =    {'febg7_u000_009';'febg7_u001_009';   'febg8_u001_008';'febg8_u002_006'};

data(3).age = 'P36-37';
data(3).expList{1,1} =   {'febg9_u000_006';                    'febh5_u000_003';                     'febl0_u000_010'                                  }; %Only pre-any cooling
data(3).expList{2,1} =   {'febg9_u000_006';'febg9_u001_000';   'febh5_u000_003';'febh5_u001_004';    'febl0_u000_010';'febl0_u001_006';'febl0_u001_018'};
data(3).expList{3,1} =   {'febg9_u000_007';'febg9_u001_004';   'febh5_u000_006';'febh5_u001_008';    'febl0_u000_012';'febl0_u001_010';'febl0_u001_019'};
data(3).expList{4,1} =   {'febg9_u000_011';'febg9_u001_008';   'febh5_u000_009';'febh5_u001_015';    'febl0_u000_015';'febl0_u001_016'                 };

data(4).age = 'Adults';
data(4).expList{1,1} =   {'febe0_u000_000';'febe0_u001_001';                    'febe2_u000_000';   'febf6_u000_000';'febf6_u001_000';                    'febf7_u000_000';                    'febg0_u000_001';'febg0_u001_000';'febg0_u002_000';'febg0_u003_002';                       'febg1_u000_000';'febg1_u001_000';    'febg2_u001_003';                                       'febg3_u000_000'}; %Only pre-any cooling
data(4).expList{2,1} =   {'febe0_u000_000';'febe0_u001_001';'febe0_u002_000';   'febe2_u000_000';   'febf6_u000_000';'febf6_u001_000';                    'febf7_u000_000';'febf7_u001_001';   'febg0_u000_001';'febg0_u001_000';'febg0_u002_000';'febg0_u003_002';                       'febg1_u000_000';'febg1_u001_000';    'febg2_u001_003';                                       'febg3_u000_000';                 'febg3_u001_011'};
data(4).expList{3,1} =   {                 'febe0_u001_004';'febe0_u002_003';   'febe2_u000_004';                    'febf6_u001_003';'febf6_u002_004';   'febf7_u000_004';                                                                       'febg0_u003_003';                                        'febg1_u001_003';    'febg2_u001_005';'febg2_u001_016';'febg2_u001_020';                      'febg3_u001_003';'febg3_u001_013'};
data(4).expList{4,1} =   {                                                                                           'febf6_u001_007';'febf6_u002_008';   'febf7_u000_007';'febf7_u001_003';                                                      'febg0_u003_004';                                        'febg1_u001_004';    'febg2_u001_013';'febg2_u001_017';'febg2_u001_021';                      'febg3_u001_004';'febg3_u001_017'};

for ageGroup = 1:4
for c = conds
    sums = [];

for e = 1:length(data(ageGroup).expList{c,1})

    exptName = data(ageGroup).expList{c,1}{e};
    animal = exptName(1:5);
    unit = exptName(8:10);
    expt = exptName(12:14);
    disp(exptName)
    load(fullfile(physDir,animal,exptName,[exptName '_id.mat']))

    linStyle = {'-','--',':'};
    if strcmp(area,'V1')
        color = 'b';
    elseif strcmp(area,'PSS')
        color = 'r';
    end


    [sumStats,~] = plotUnits(animal,unit,expt,area,anaMode,'ranksum',0.01,0,0,0,dataFold);
    close all
    if isempty(sumStats)
        continue
    end
    sumStats.dsi(sumStats.dsi>1) = 1;
    sumStats = sumStats(sumStats.goodUnit==1,:);
    sums = vertcat(sums,sumStats);

end

data(ageGroup).sums{c,1} = sums;


end


end

%% Plotting


figure;
for a = 2:4
    subplot(1,3,a-1);hold on
    clear x y
for c = conds

    x = mean(vertcat(data(a).sums{c}.tuningX{:}));
%     y = vertcat(data(a).sums{c}.tuningY{:});
    for u = 1:height(data(a).sums{c})

        y(u,:) = data(a).sums{c}.tuningY{u};
        y(u,:) = y(u,:)/max(y(u,:));

    end

    sem = std(y)./sqrt(size(y,1));
    plot(x,mean(y),condColor{c},'LineWidth',2);
    plot(repmat(x,2,1),mean(y)+([-1;1]*sem),condColor{c},'LineWidth',2)
    title(data(a).age)

end
xticks([-180:90:180])
ylabel('normalized firing rate')
end




figure;
for a = 2:4
    subplot(1,3,a-1);hold on
    clear x y
for c = conds

    x = mean(vertcat(data(a).sums{c}.tuningX{:}));
    y = vertcat(data(a).sums{c}.tuningY{:});
    

    sem = std(y)./sqrt(size(y,1));
    plot(x,mean(y),condColor{c},'LineWidth',2);
    plot(repmat(x,2,1),mean(y)+([-1;1]*sem),condColor{c},'LineWidth',2)
    title(data(a).age)

end
xticks([-180:90:180])
ylabel('firing rate')
end




figure;hold on

tiledlayout(3,4)
for t = 1:3
    
    for a = 1:length(data)
    
    nexttile;hold on
    for c = conds
        
        if t == 1
            dist{t,a,c} = data(a).sums{c,1}.rPref;
            yLbl = 'rPref';
            lims = [1 100 1000];
            offsetScale = 3;
            dpScale = 20;
        elseif t == 2
            dist{t,a,c} = data(a).sums{c,1}.dsi;
            yLbl = 'DSI';
            lims = [0 1 1000];
            offsetScale = 3;
            dpScale = 20;
        elseif t == 3
            dist{t,a,c} = 1-data(a).sums{c,1}.dcv;
            yLbl = '1-DCV';
            lims = [0 1 1000];
            offsetScale = 3;
            dpScale = 20;
        end

        cdf = cdfplot(dist{t,a,c});
        cdf.Color = condColor{c};
        ylabel('percentile')
        xlabel(yLbl)

%         histogram(dist{t,a,c},'FaceColor',condColor{c},'FaceAlpha',0.5,'EdgeColor','none')


%         expts = unique(data(a).sums{c,1}.exptName);
%         nE = length(expts);
%         for e = 1:nE
%             binID = strcmp(data(a).sums{c,1}.exptName,expts{e});
%             exptID(binID) = e;
%         end
% 
%         if c == 1
%             
%         elseif c == 2
% 
%         end
% 
%         rng = lims(2)-lims(1);
%         linWdth = rng/(1*10^1);
%         xBar = mean(dist{t,a,c});
%         med = median(dist{t,a,c});
%         sem = std(dist{t,a,c})/sqrt(length(dist{t,a,c}));
%         v = violin(dist{t,a,c},0.4,lims,1,0);
%         z = c*offsetScale;
%         s = offsetScale/dpScale;
% %         plot(((exptID-(nE/2)-0.5)*s*indAnimals)+z,dist{t,a,c}','.','Color',condColor{c});
%     %     plot((vertcat(zeros(size(dist')),ones(size(dist')))*(s))+(z),repmat(dist',2,1),'k');
%     %     for dp = 1:length(dist)
%     %         plot(((v.g(dp,:)*-1)*(s))+(z),v.gX(1,:),'b')
%     %     end
%         patch([v.kdeX -1*fliplr(v.kdeX)]+z,[v.kdeY fliplr(v.kdeY)],condColor{c},'FaceAlpha',0.1,'EdgeColor','none')
%         plot(v.kdeX+z,v.kdeY,'Color',condColor{c})
%         plot((v.kdeX*-1)+z,v.kdeY,'Color',condColor{c})
%         plot([z+s z-s],repmat(xBar,2,1),'r','lineWidth',2)
%         plot([z+s z-s],repmat(med,2,1),'b','lineWidth',2)
%     %     plot([z z],xBar+([-1 1]*sem),'r')
%         
%         clear exptID binID


    end

%     ylim([lims(1) lims(2)])
%     ylabel(yLbl)
%     xlim([offsetScale 5*offsetScale])
%     xticks(conds*offsetScale)
%     xticklabels({'before','during','after'})
% %     set(gcf, 'Position',  [0, 100, 300, 300])
    
    end

end









