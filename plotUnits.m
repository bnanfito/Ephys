%plotUnits3

function [sumStats,spks] = plotUnits(animal,unit,expt,probe,anaMode,visTest,alpha,plt,saveSum,saveFigs,dataFold)

% close all
% clear all
% animal = 'FEAO4';
% unit = '001';
% expt = '000';
% probe = 1;
% anaMode = 'MU';
% onlyGoodUnits = 0;
% saveSum = 0;
% saveFigs = 0;
% plt = 0;

split = 0; splitInd = 35;
polar = 0;
alignBit = 1;
sortTrialBit = 0;
stat = 'sem';
baseName = [animal '_u' unit '_' expt];

physDir = fullfile(dataFold,'Ephys',animal,baseName); 
figDir = fullfile(dataFold,'Figures',animal,baseName); 
sumDir = fullfile(dataFold,'SummaryStats',animal,baseName); 

load(fullfile(physDir,[baseName '_id.mat']),'id')
load(fullfile(physDir,[baseName '_trialInfo']),'trialInfo')
load(fullfile(physDir,[baseName '.analyzer']),'-mat')

varInfo = whos('probe');
if strcmp(varInfo.class,'char')
    area = probe;
    probe = find(strcmp({id.probes(:).area}',area));
end

area = id.probes(probe).area;
sf = id.sampleFreq;
% kernel = ones(1,0.1*sf)*(1/(0.1*sf));
kernel = normpdf(-3:6/2000:3);
predelay = getparam('predelay',Analyzer);
stimTime = getparam('stim_time',Analyzer);
postdelay = getparam('postdelay',Analyzer);
trialL = predelay+stimTime+postdelay;
nTrials = length(trialInfo.triallist);
nConds = length(unique(trialInfo.triallist));
nDom = length(trialInfo.dom);
nReps = nTrials/nConds;
if isempty(trialInfo.blankId)
    blank = zeros(1,nConds)==1;
else
    blank = (1:nConds)==trialInfo.blankId;
end
[~,sortTrialInd] = sort(trialInfo.triallist);
nT = ((1:(trialL*sf))-(predelay*sf))/sf;
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};

for c = 1:nConds

    trialRepMap(:,c) = find(trialInfo.triallist == c);

end
clear c


%% DATA PROCESSING
% collect spike data from either spkSort file (if single unit analysis) or
% MUspkMerge file (if mua). Determine which units pass inclusion criteria
% (are 'goodUnits')


[spks,trialExclude] = orgSpks(animal,unit,expt,probe,anaMode,dataFold);
[goodUnits,pVis] = screenUnits(spks,anaMode,blank,visTest,alpha);


if sum(goodUnits) == 0
disp('no good units')
end
uID = vertcat(spks.unitId);

if split == 1
    half{1} = 1:splitInd-1;
    half{2} = splitInd:nTrials;
end


%% PLOT
sdf = zeros(size(spks(1).train,1),length(spks));
for u = 1:length(spks)
    
    exptNames{u,1} = baseName;
    probeIDs(u,1) = probe;
    areaIDs{u,1} = area;

    if plt==1

    figure;hold on
    % RASTER PLOT
    subplot(1,3,1);hold on
%     [x,y] = find(spks(u).train);
    x = spks(u).stimCent(1,:);
    y = spks(u).stimCent(2,:);
    if sortTrialBit == 1
        [y,~] = find(sortTrialInd==y');
    end
    patch([0 stimTime stimTime 0],[0 0 nTrials+1 nTrials+1],'k','EdgeColor','none','FaceAlpha',0.2)
    for t = 1:nTrials
        if ismember(t,find(trialExclude))
            if sortTrialBit == 1
                [t,~] = find(sortTrialInd==t);
            end
            patch([-predelay stimTime+postdelay stimTime+postdelay -predelay],[t-0.5 t-0.5 t+0.5 t+0.5],'r','EdgeColor','none','FaceAlpha',0.2)
        end
    end
    plot(x,y,'k.')
    axis tight
    xlim([-predelay stimTime+postdelay])
    xlabel('time (s) relative to stim onset')
    ylabel('trial #')
    
    if split == 1
    for h = 1:2
        plot(x(ismember(y,half{h})),y(ismember(y,half{h})),'o','MarkerSize',5,'Color',colors{h})
    end
    end
    
    % SDF
    subplot(1,3,2);hold on
    sdf(:,u) = conv(sum(spks(u).train,2),kernel,'same');
    patch([0 stimTime stimTime 0],[0 0 max(sdf(:,u))+(max(sdf(:,u))*0.1) max(sdf(:,u))+(max(sdf(:,u))*0.1)],'k','EdgeColor','none','FaceAlpha',0.2)
    plot(nT,sdf(:,u),'k')
    if split == 1
    for h = 1:2
        sdfH{h}(:,u) = conv(sum(spks(u).train(:,half{h}),2),kernel,'same');
        plot(nT,sdfH{h}(:,u),'--','Color',colors{h})
    end
    end
    axis tight
    xlabel('time (s) relative to stim onset')
    ylabel('SDF')
    end
    
    % TUNING CURVE
    if (length(trialInfo.dom)==1 && strcmp(trialInfo.dom{1},'ori')) % for experiments where only orientation/direction of motion is changing


        c = trialInfo.domval; % vectors of conditions (x axis of tuning curve; e.g. orientations)
        r = spks(u).fr.bc(:,~blank);
        rBlank = spks(u).fr.stim(:,blank);
        rB(u,1) = mean(rBlank);
        sem = std(r,'omitnan')/sqrt(size(r,1)); % nConds long vector of standard error of the mean for each condition
        ci95 = confInt(r); % nConds long vector of 95% confidence interval for each condition
        
        if split == 1
        for h = 1:2
            if h==1
                hInd = trialRepMap<splitInd;
            elseif h==2
                hInd = trialRepMap>=splitInd;
            end
            rH = r; rH(~hInd)=nan;
            yH{h} = mean(rH,'omitnan');
            rPrefH{h}(u,1) = max(yH{h});
        end
        end

        rPref = max(mean(r,'omitnan')); rP(u,1) = rPref;
        cPref = find(mean(r,'omitnan')==rPref);

        if length(c) == 2 %%% bidirectional training %%%

            if length(cPref)>1 % if there is more than one pk with rPref
                cPref = cPref(1);
            end
            cP(u,1) = c(cPref);
            
            x = c;                  tuningX{u,1} = x;
            y = mean(r,'omitnan');  tuningY{u,1} = y; 

            dpi(u,1) = abs(diff(y))/max(y);

            if  split == 1
            for h = 1:2
                dpiH{h}(u,1) = abs(diff(yH{h}))/max(yH{h});
            end
            end

            if plt == 1
                subplot(1,3,3);hold on
                for o = 1:2
                    plot(c(o)+(randn(size(r,1),1)*10),r(:,o),'.','Color',[0.8 0.8 0.8])
                    plot(repmat(c(o),2,1),y(o)+(sem(o)*[-1;1]),'r','LineWidth',2)
                end
                plot(x,y,'-o','Color','k','LineWidth',2)

                if split == 1
                for h = 1:2
                    plot(x,yH{h},'--o','Color',colors{h},'LineWidth',2)
                end
                end
    
                xlim([min(c)-50 max(c)+50])
                xticks(c)
                xlabel('direction of motion (deg)')
                ylabel('BCFR (Hz.)')
                title(['dpi = ' num2str(dpi(u,1))])
            end


        elseif length(c) > 2 %%% ori12 / ori16 experiment %%%

            % calculare summary metrics in ori space (mod 180)
            c_ori = mod(c,180);
            oris = unique(c_ori);
            for o = 1:length(oris)
                rTemp = r(:,c_ori == oris(o));
                r_ori(:,o) = rTemp(:);
                clear rTemp
            end
            rPref_ori = max(mean(r_ori,'omitnan'));
            cPref_ori = oris(mean(r_ori,'omitnan')==rPref_ori);
            if length(cPref_ori)>1 % if there is more than one pk with rPref
                rIn = mean(r_ori,'omitnan');
                pks = rIn == rPref_ori;
                rConv = rIn([end 1:end 1]);
                rConv = conv(rConv,ones(1,3)*(1/3),'same');
                rConv = rConv(1+1:end-1);
                rConv(~pks) = 0;
                cPref_ori = oris(find(rConv==max(rConv),1,'first'));
                clear rIn pks rConv 
            end
            cNull_ori = mod(cPref_ori+90,180);
            rNull_ori = mean(r_ori(:,oris==cNull_ori),'omitnan');
            osi(u,1) = abs(rPref_ori-rNull_ori)/rPref_ori;

            % calculate summary metrics in dir space (mod 360)
            if length(cPref)>1 % if there is more than one pk with rPref
                rIn = mean(r,'omitnan');
                pks = rIn == rPref;
                rConv = rIn([end 1:end 1]);
                rConv = conv(rConv,ones(1,3)*(1/3),'same');
                rConv = rConv(1+1:end-1);
                rConv(~pks) = 0;
                cPref = find(rConv==max(rConv),1,'first');
                clear rIn pks rConv 
            end
            cP(u,1) = c(cPref);
            cNull = c==mod(cP(u,1) + 180,360); cN(u,1) = c(cNull);
            rNull = mean(r(:,c==cN(u,1)),'omitnan'); rN(u,1) = rNull;
            dsi(u,1) = abs(rPref-rNull)/rPref;
            mv = meanvec(c,mean(r,'omitnan'));
            ldir(u,1) = mv.ldir;
            lori(u,1) = mv.lori;
            
            % align tuning curves relative to cPref
            if alignBit == 1
                [x,y,cMap] = alignDirTuning(c',mean(r,'omitnan')); 
                tuningX{u,1} = x; 
                tuningY{u,1} = y;
            elseif alignBit == 0
                x = c'; y = mean(r,'omitnan');
                tuningX{u,1} = x;
                tuningY{u,1} = y;
                cMap = 1:length(c);
            end
            sem = sem(cMap);
            ci95 = ci95(:,cMap);
            if split == 1
                yH{1} = yH{1}(cMap);
                yH{2} = yH{2}(cMap);
            end

            if plt ==1
                if polar == 1
                    subplot(1,3,3,polaraxes);hold on
                else
                    subplot(1,3,3);hold on
                end
                
                patch([min(x) max(x) max(x) min(x)],[min(rBlank) min(rBlank) max(rBlank) max(rBlank)],'r','EdgeColor','none','FaceAlpha',0.1)
                yline(mean(rBlank),'r--')
                plot(x,r(:,cMap)','--.','Color',[0.8 0.8 0.8],'MarkerSize',10)
                plot(x,y,'-o','Color','k','LineWidth',2) % tuning curve (mean across reps)
                plot(repmat(x,2,1),y+([-1;1]*sem),'k','LineWidth',2)
                yline(0,'k--')
                if alignBit == 1
                    xticks([-180 -90 0 90 180])
                elseif alignBit == 0
                    xticks([0 90 180 270])
                end
                xlabel('direction of motion (deg)')
                xlim([min(x) max(x)])
                ylabel('BCFR (Hz.)')
                title(['dsi = ' num2str(dsi(u,1)) ', 1-dcv = ' num2str(1-dcv(u,1))])
            end

        end

    elseif sum(strcmp(trialInfo.dom,'y_size'))>0 % for experiments where stimulus alternates between full field and hemifield 

        disp('still editing this section')
%         oriInd = find(strcmp(trialInfo.dom,'ori'));
%         oris = unique(trialInfo.domval(:,oriInd));
%         ySizeInd = find(strcmp(trialInfo.dom,'y_size'));
%         ySizes = unique(trialInfo.domval(:,ySizeInd));
%         ff_cndInd = trialInfo.domval(:,ySizeInd) == ySizes(ySizes>180);
%         ff_trialInd = ismember(trialInfo.triallist,find(ff_cndInd));
%         hemi_cndInd = trialInfo.domval(:,ySizeInd) == ySizes(ySizes<=75);
%         hemi_trialInd = ismember(trialInfo.triallist,find(hemi_cndInd));
%         rBlank = spks(u).fr.bc(:,blank);
%         rB(u,1) = mean(rBlank,'omitnan');
% 
%         for stim = 1:2
% 
%             if stim == 1
%                 idxCnd = ff_cndInd;
%                 idxTrial = ff_trialInd;
%                 clr = 'g';
%             elseif stim == 2
%                 idxCnd = hemi_cndInd;
%                 idxTrial = hemi_trialInd;
%                 clr = 'k';
%             end
% 
%             c{stim} = trialInfo.domval(idxCnd,oriInd);
%             r{stim} = spks(u).fr.bc(:,idxCnd);
%             rPref{stim} = max( mean(r{stim}) );
%             cPref{stim} = find(mean(r{stim},'omitnan') == rPref{stim});
%             if length(cPref{stim})>1 % if there is more than one pk with rPref
%                 rIn = mean(r{stim},'omitnan');
%                 pks = rIn == rPref{stim};
%                 rConv = rIn([end 1:end 1]);
%                 rConv = conv(rConv,ones(1,3)*(1/3),'same');
%                 rConv = rConv(1+1:end-1);
%                 rConv(~pks) = 0;
%                 cPref{stim} = find(rConv==max(rConv),1,'first');
%                 clear rIn pks rConv 
%             end
%             cNull{stim} = mod(cPref{stim}+180,360);
%             rNull{stim} = r{stim}(c{stim}==cNull{stim});
%             dsi_tmp{stim} = abs(rPref{stim}-rNull{stim})/rPref{stim};
%             mv = meanvec(c{stim},mean(r{stim},'omitnan'));
%             ldir_tmp{stim} = mv.ldir;
%             lori_tmp{stim} = mv.lori;
%             
% 
% 
% 
%             if plt == 1
%                 subplot(1,3,3);hold on
%                 plot(c{stim},mean(r{stim},'omitnan'),'-o','LineWidth',2,'Color',clr)
%                 plot(c{stim},r{stim}','--','LineWidth',1,'Color',clr)
% 
%             end
%             
% 
% 
%         end
% 
%         cP{u,1} = c{1}([cPref{:}]);
%         rP{u,1} = [rPref{:}];
%         rN{u,1} = [rNull{:}];
%         cN{u,1} = c{1}([cPref{:}]);
%         dsi{u,1} = [dsi_tmp{:}];
%         ldir{u,1} = [ldir_tmp{:}];
%         lori{u,1} = [lori_tmp{:}];
%         tuningX{u,1} = c;
%         tuningY{u,1} = r;



    elseif length(trialInfo.dom) == 2  % for experiments where another variable changes in addition to ori/dir (e.g. contrast experiments, spatial frequency, etc)

        oriInd = find(strcmp(trialInfo.dom,'ori'));
        notOriInd = find(~((1:2)==oriInd));
        oris = unique(trialInfo.domval(:,oriInd));
        rBlank = spks(u).fr.bc(:,blank);
        rB(u,1) = mean(rBlank);
        for o = 1:length(oris)
            oriIdx = trialInfo.domval(:,oriInd)==oris(o);
            c{o} = trialInfo.domval(oriIdx,notOriInd);
            r{o} = mean(spks(u).fr.bc(:,oriIdx),'omitnan');
            oriMax(o) = max(r{o});
%             [nakaX,nakaY,~,~,~,~,~] = nakaRush(r{o},c{o},1,0);
            if plt == 1
                subplot(1,3,3);hold on
                plot(c{o},r{o},'-o','LineWidth',2)
%                 plot(nakaX,nakaY,'--','LineWidth',2)
            end

        end
        if length(unique(oriMax)) == 1   
            oriPref = 1;
        elseif length(unique(oriMax)) > 1
            oriPref = find(oriMax == max(oriMax));
        end
        tuningX{u,1} = c{oriPref};
        tuningY{u,1} = r{oriPref};
        cPref = find(trialInfo.domval(:,oriInd) == oris(oriPref) & trialInfo.domval(:,notOriInd) == 100);
        rP(u,1) = oriMax(oriPref);
        cP(u,1) = c{oriPref}( find(r{oriPref} == rP(u,1),1) );

    end

    %Calculate latency
    binWidth=0.010; %sec
    % 1sec baseline and 1sec stim period
    startBin=ceil(-1/binWidth)*binWidth; %need multiple of binWidth to make 0 an edge
    stopBin=floor(1/binWidth)*binWidth;
    binVec=[startBin:binWidth:stopBin];
    if length(cPref)==1
        prefTrials = find(trialInfo.triallist==cPref);
    else
        prefTrials = find(trialInfo.triallist==cPref{1});
    end
    for rep = 1:length(prefTrials)
        t = prefTrials(rep);
        spkTs{rep} = spks(u).stimCent(1,spks(u).stimCent(2,:)==t);
        N(rep,:) = histcounts(spkTs{rep},binVec);
    end
    avgN = mean(N,1)/binWidth;
    avgBase=mean(avgN(binVec<0));
    stdBase=std(avgN(binVec<0));
    cumN=cumsum(avgN-avgBase); %cumsum(1): value 1 in input
    diffN=diff(cumN);
    %criterion 1: cumsum over threshold and 2 increasing bins
    idx=find(cumN(1:end-2)>2*stdBase & diffN(1:end-1)>0 & diffN(2:end)>0 ...
        & binVec(1:end-3)>0,1);
    if ~isempty(idx)
        latCh(u,1)=binVec(idx);
    else
        latCh(u,1)=NaN;
    end

    %criterion 2: cumsum over threshold and 3 increasing bins
    idx=find(cumN(1:end-3)>2*stdBase & ...
        diffN(1:end-2)>0 & diffN(2:end-1)>0 & diffN(3:end)>0 ...
        & binVec(1:end-4)>0,1);
    if ~isempty(idx)
        latCh2(u,1)=binVec(idx);
    else
        latCh2(u,1)=NaN;
    end



     if plt == 1
        ttl = [baseName ' p' num2str(probe) ' (' id.probes(probe).area ') ' anaMode '#' num2str(uID(u))];
        if ~ismember(u,find(goodUnits))
            ttl = [ttl ' (BAD UNIT)'];
        end
        sgtitle(ttl)
        set(gcf, 'Position',  [0, 100, 1000, 300])
    
        if saveFigs == 1
            if ~isdir(figDir)
                mkdir(figDir)
            end
            saveas(gcf,fullfile(figDir,[baseName '_p' num2str(probe) '_' anaMode num2str(u)]),'fig')
            saveas(gcf,fullfile(figDir,[baseName '_p' num2str(probe) '_' anaMode num2str(u)]),'svg')
            saveas(gcf,fullfile(figDir,[baseName '_p' num2str(probe) '_' anaMode num2str(u)]),'jpeg')
        end
    
     end

end

%% Compute F1/F0


%% Make Summary Stats Table

varNames = {'exptName','probe','area','uID','goodUnit','pVis','fr','rPref','cPref','rBlank','lat1','lat2'};
sumStats = table(exptNames,probeIDs,areaIDs,uID,goodUnits',pVis,vertcat(spks.fr),rP,cP,rB,latCh,latCh2,'VariableNames',varNames);
sumStats.tuningX = tuningX;
sumStats.tuningY = tuningY;
if strcmp(anaMode,'MU')
    sumStats.xPos = vertcat(spks.xPos);
    sumStats.zPos = vertcat(spks.zPos);
end
if exist('rN','var')
    sumStats.rNull = rN;
    sumStats.cNull = cN;
end
if exist('dsi','var')
    sumStats.dsi = dsi;
    sumStats.ldir = ldir;
    sumStats.osi = osi;
    sumStats.lori = lori;
end
if exist('dpi','var')
    sumStats.dpi = dpi;
end
if exist('c50','var')
    sumStats.c50 = c50;
end

if saveSum == 1 
    if ~isfolder(sumDir)
        mkdir(sumDir)
    end
    save(fullfile(sumDir,[baseName '_p' num2str(probe) '_sumStats' anaMode '.mat']),'sumStats','spks')
end


%% Plot Summary

if plt == 1 
    goodInd = sumStats.goodUnit == 1;
    if exist('dsi','var')
        
        figure
    
        subplot(2,3,1); hold on
        if sum(goodInd)>0
            sumStats.rPref(goodInd)
            p1 = cdfplot(sumStats.rPref(goodInd)); p1.Color = 'k'; p1.LineWidth = 2;
        end
        p2 = cdfplot(sumStats.rPref); p2.LineStyle = '--'; p2.Color = 'k'; p2.LineWidth = 2;
        xlabel('response to pref. (BCFR in HZ)');
        ylabel('percentile');
        legend({['good units (n=' num2str(length(find(goodUnits))) ')'],'all units'})
    
        subplot(2,3,2); hold on
        d = sumStats.dsi(goodInd); d(d>1) = 1;
        if sum(goodInd)>0
            p3 = cdfplot(d); p3.Color = 'k'; p3.LineWidth = 2;
        end
        p4 = cdfplot(sumStats.dsi); p4.LineStyle = '--'; p4.Color = 'k'; p4.LineWidth = 2;
        xlabel('dsi'); xlim([0 1])
        ylabel('percentile')
    
        subplot(2,3,3); hold on
        if sum(goodInd)>0
            p5 = cdfplot(1-sumStats.dcv(goodInd)); p5.Color = 'k'; p5.LineWidth = 2;
        end
        p6 = cdfplot(1-sumStats.dcv); p6.LineStyle = '--'; p6.Color = 'k'; p6.LineWidth = 2;
        xlabel('1-dcv'); xlim([0 1])
        ylabel('percentile')

        subplot(2,3,4); hold on
        x = mean(vertcat(sumStats.tuningX{:}));
        y = vertcat(sumStats.tuningY{:});
        sem = std(y)/sqrt(size(y,1));
        plot(x,mean(y),'k--o','LineWidth',2)
        plot(repmat(x,2,1),mean(y)+([-1;1]*sem),'k','LineWidth',2)
        y = vertcat(sumStats(goodInd,:).tuningY{:});
        if size(y,1)==1
            sem = std(y)/sqrt(1);
            plot(x,y,'k-o','LineWidth',2)
            plot(repmat(x,2,1),y+([-1;1]*sem),'k','LineWidth',2)
        else
            sem = std(y)/sqrt(size(y,1));
            plot(x,mean(y),'k-o','LineWidth',2)
            plot(repmat(x,2,1),mean(y)+([-1;1]*sem),'k','LineWidth',2)
        end
        if alignBit == 1
            xticks([-180 -90 0 90 180])
        else
            xticks([0 90 180 270])
        end
    
        subplot(2,3,5); hold on
        x = sumStats.dsi(goodInd); x(x>1) = 1;
        y = sumStats.rPref(goodInd);
        id = sumStats.uID(goodInd);
        s1 = plot(x,y,'k.','MarkerSize',20);
        if ~isempty(y)
            s1.DataTipTemplate.DataTipRows(1:3) = [dataTipTextRow('dsi = ',x) dataTipTextRow('rPref = ',y) dataTipTextRow([anaMode ' = '],id)];
        end
        x = sumStats.dsi(~goodInd); x(x>1) = 1;
        y = sumStats.rPref(~goodInd);
        id = sumStats.uID(~goodInd);
        s2 = plot(x,y,'ko');
        if ~isempty(y)
            s2.DataTipTemplate.DataTipRows(1:3) = [dataTipTextRow('dsi = ',x) dataTipTextRow('rPref = ',y) dataTipTextRow([anaMode ' = '],id)];
        end
        legend({['good units (n=' num2str(length(find(goodUnits))) ')'],'bad units'})
        xlabel('dsi'); xlim([0 1])
        ylabel('rPref')
        
    
        subplot(2,3,6); hold on
        x = 1-sumStats.dcv(goodInd);
        y = sumStats.rPref(goodInd);
        id = sumStats.uID(goodInd);
        s3 = plot(x,y,'k.','MarkerSize',20);
        if ~isempty(y)
            s3.DataTipTemplate.DataTipRows(1:3) = [dataTipTextRow('1-dcv = ',x) dataTipTextRow('rPref = ',y) dataTipTextRow([anaMode ' = '],id)];
        end
        x = 1-sumStats.dcv(~goodInd);
        y = sumStats.rPref(~goodInd);
        id = sumStats.uID(~goodInd);
        s4 = plot(x,y,'ko');
        if ~isempty(y)
            s4.DataTipTemplate.DataTipRows(1:3) = [dataTipTextRow('1-dcv = ',x) dataTipTextRow('rPref = ',y) dataTipTextRow([anaMode ' = '],id)];
        end
        legend({['good units (n=' num2str(length(find(goodUnits))) ')'],'bad units'})
        xlabel('1-dcv'); xlim([0 1])
        ylabel('rPref')
    
    end
    
    if exist('dpi','var')
        
        figure
    
        subplot(2,2,1); hold on
        if sum(goodInd)>0
            p1 = cdfplot(sumStats.rPref(sumStats.goodUnit)); p1.Color = 'k'; p1.LineWidth = 2;
        end
        p2 = cdfplot(sumStats.rPref); p2.LineStyle = '--'; p2.Color = 'k'; p2.LineWidth = 2;
        if split == 1
        for h = 1:2
            pH = cdfplot(rPrefH{h}(sumStats.goodUnit)); pH.Color = colors{h}; pH.LineWidth = 2;
        end
        end
        xlabel('response to pref. (BCFR in HZ)');
        ylabel('percentile');
        legend({['good units (n=' num2str(length(find(goodUnits))) ')'],'all units'})
    
        subplot(2,2,2); hold on
        if sum(goodInd)>0
            p3 = cdfplot(sumStats.dpi(sumStats.goodUnit==1)); p3.Color = 'k'; p3.LineWidth = 2;
        end
        p4 = cdfplot(sumStats.dpi); p4.LineStyle = '--'; p4.Color = 'k'; p4.LineWidth = 2;
        if split == 1
        for h = 1:2
            pH = cdfplot(dpiH{h}(sumStats.goodUnit)); pH.Color = colors{h}; pH.LineWidth = 2;
        end
        end
        xlabel('dpi')
        xlim([0 1])
        ylabel('percentile')
    
        subplot(2,2,4); hold on
        x = sumStats.dpi(goodInd);
        y = sumStats.rPref(goodInd);
        id = sumStats.uID(goodInd);
        s1 = plot(x,y,'k.','MarkerSize',20);
        if ~isempty(y)
            s1.DataTipTemplate.DataTipRows(1:3) = [dataTipTextRow('dpi = ',x) dataTipTextRow('rPref = ',y) dataTipTextRow([anaMode ' = '],id)];
        end
        x = sumStats.dpi(~goodInd);
        y = sumStats.rPref(~goodInd);
        id = sumStats.uID(~goodInd);
        s2 = plot(x,y,'ko');
        if ~isempty(y)
            s2.DataTipTemplate.DataTipRows(1:3) = [dataTipTextRow('dpi = ',x) dataTipTextRow('rPref = ',y) dataTipTextRow([anaMode ' = '],id)];
        end
        xlabel('dpi')
        xlim([0 1])
        ylabel('repsonse to pref. (BCFR in Hz.)')
        legend({['good units (n=' num2str(length(find(goodUnits))) ')'],'bad units'})
        
    
    end
    
    if exist('c50','var')
    
    
    end
    
    
%     set(gcf,'Position',[0 0 1000 1000])
    sgtitle([animal ' u' unit ' e' expt ' p' num2str(probe) '(' area ') summary plot'])
    
    if saveFigs == 1
        if ~isfolder(figDir)
            mkdir(figDir)
        end
        saveas(gcf,fullfile(figDir,[baseName '_p' num2str(probe) '_' anaMode 'sumPlot' ]),'fig')
        saveas(gcf,fullfile(figDir,[baseName '_p' num2str(probe) '_' anaMode 'sumPlot']),'svg')
        saveas(gcf,fullfile(figDir,[baseName '_p' num2str(probe) '_' anaMode 'sumPlot']),'jpeg')
    end
end

end