%Bayesian Decoder
%written by Brandon Nanfito

clear all
close all

anaMode = 'SU';
trainProj = 'Train_V1Cool';
area = 'PSS';
trainC = 'af';
% dataFold = 'F:\Brandon\data';
dataFold = '/Volumes/NielsenHome2/Brandon/data';
dsFold = fullfile(dataFold,'dataSets','training',trainProj,anaMode);
load(fullfile(dsFold,[trainProj '_' anaMode 'dataSet.mat']))

if strcmp(area,'V1')
    if strcmp(trainC,'bf')
        tbl = data.v1bf;
    elseif strcmp(trainC,'af')
        tbl = data.v1af;
    end
elseif strcmp(area,'PSS')
    if strcmp(trainC,'bf')
        tbl = data.pssbf;
    elseif strcmp(trainC,'af')
        tbl = data.pssaf;
    end
end
tbl = tbl(tbl.goodUnit,:);
nU = height(tbl);

expts = unique(tbl.exptName);
for e = 1:length(expts)
    eIdx = find(strcmp(tbl.exptName,expts{e}));
    eIdxStart(e) = eIdx(1);
    eIdxEnd(e) = eIdx(end);
end

[~,sortIdx] = sort(tbl.oriPref);
% tbl = tbl(sortIdx,:);

R = cat(3,tbl.response{:});
nReps = size(R,1);
nConds = size(R,2);
Rt = reshape(R,nReps*nConds,nU);
Rt(Rt<0) = 0;
Rtnorm = Rt./max(Rt);

figure;

subplot(2,2,1)
imagesc(Rt)
xticks(eIdxEnd)
title('Response Matrix')
xlabel('unit')
ylabel('trial')

subplot(2,2,2)
imagesc(Rtnorm)
xticks(eIdxEnd)
title('norm. Response Matrix')
xlabel('unit')
ylabel('trial')

subplot(2,2,3)
imagesc(Rtnorm(:,sortIdx))
title('norm. Response Matrix')
xlabel('unit (sort by oriPref)')
ylabel('trial')

sgtitle([trainProj ' ' anaMode ' ' area ' ' trainC])
figName = [trainProj '_RM_' anaMode '_' area '_' trainC];
saveas(gcf,fullfile(dsFold,figName),'fig')







% figure;hold on
% for i = 1:4
% 
%     if i==1
%         I = 'v1bf';
%         tbl = data.v1bf;
%         clr = 'b';
%         linStyl = '--';
%     elseif i==2
%         I = 'v1af';
%         tbl = data.v1af;
%         clr = 'b';
%         linStyl = '-';
%     elseif i==3
%         I = 'pssbf';
%         tbl = data.pssbf;
%         clr = 'r';
%         linStyl = '--';
%     elseif i==4
%         I = 'pssaf';
%         tbl = data.pssaf;
%         clr = 'r';
%         linStyl = '-';
%     end
%     
%     nU = height(tbl);
%     for u = 1:nU
%     
%         ldr{i}(u) = tbl.ldr(u);
%         R = tbl.response{u};
%         rVar{i}(u) = mean(var(R,'omitnan'));
%         nStim = size(R,2);
%         nRep = size(R,1);
%         nTrial = nStim*nRep;
%         R = R(:);
%         meanR{i}(u) = mean(R,'omitnan');
%         if strcmp(anaMode,'MU')
%             skipTrial = isnan(R);
%         end
%         S = repmat(1:nStim,nRep,1);
%         S = S(:);
%         for t = 1:nTrial
%     
%             Pr = sum(R==R(t))/nTrial;
%             for s = 1:nStim
%                 Rs = R(S==s);
% 
%                 Ps = 1/nStim;
% %                 Pr_s = sum(Rs==R(t))/nRep;
%                 Pr_s = gauss(R(t),mean(Rs),std(Rs));
%                 Gs{i,s,u} = vertcat( gauss(linspace(min(R),max(R)),mean(Rs),std(Rs)), ...
%                                       linspace(min(R),max(R)) );
%     
%                 Ps_r{i}(s,t,u) = (Pr_s*Ps)/(Pr);
%             end
%     
%             guess = find(Ps_r{i}(:,t,u)==max(Ps_r{i}(:,t,u)));
%             if length(guess)>1
%                 guess = guess(randi(length(guess),1));
%             end
%             
%             if (strcmp(anaMode,'MU') && ismember(t,find(skipTrial))) || isempty(guess)
%                 sP{i}(t,u) = nan;
%             else
%                 sP{i}(t,u) = guess;
%             end
%     
%         end
%     
% %         if strcmp(anaMode,'MU')
% %             S = S(~skipTrial);
% %             sP{i}(:,u) = sP{i}(~skipTrial,u);
% %         end
% 
%         %confusion matrix
%         for sA = 1:nStim
%             for sB = 1:nStim
% 
%                 cMat{i}(sA,sB,u) = sum((S==sA) & (sP{i}(:,u)==sB));
% 
%             end
%         end
% 
%         hits = S==sP{i}(:,u);
%         hitRate{i}(u) = sum(hits)/nTrial;
% 
%         %calculate p-value
%         pVal{i}(u) = 0;
%         k = sum(hits); %number of hits
%         n = length(hits); %number of trials
%         p = Ps; %probability of stimuli
%         for j = k:n
%             Pj = nchoosek(n,j) *  (p^j) * ((1-p)^(n-j));
%             pVal{i}(u) = sum([pVal{i}(u),Pj]);
%         end
% 
% %         if ismember(u,find(tbl.goodUnit))
% %         figure
% %             subplot(2,2,1)
% %             histogram(R)
% %             xline(meanR{i}(u),'r--')
% %             ylabel('count')
% %             xlabel('spk count')
% %     
% %             for s = 1:nStim
% %                 subplot(nStim,2,2*s);hold on
% %                 histogram(R(S==s),min(R)-1:1:max(R)+1)
% %                 plot(Gs{i,s,u}(2,:),Gs{i,s,u}(1,:),'r-')
% %                 ylabel(num2str(s))
% %                 xlim([min(R)-1 max(R)+1])
% %             end
% % 
% %             subplot(2,2,3)
% %             imagesc(cMat{i}(:,:,u))
% %             ylabel('actual stim.')
% %             xlabel('predicted stim.')
% %     
% %             sgtitle([I ' unit#' num2str(u) '; perf = ' num2str(hitRate{i}(u)) '; pVal = ' num2str(pVal{i}(u))])
% %             fName = fullfile(dsFold,['bayesDecoder_' I anaMode num2str(u) '.fig']);
% %             saveas(gcf,fName)
% %         end
% 
%     
%     end
%     
%     hitRate{i} = hitRate{i}(tbl.goodUnit);
%     rVar{i} = rVar{i}(tbl.goodUnit);
%     pVal{i} = pVal{i}(tbl.goodUnit);
%     cMat{i} = cMat{i}(:,:,tbl.goodUnit);
%     meanR{i} = meanR{i}(tbl.goodUnit);
%     ldr{i} = ldr{i}(tbl.goodUnit);
%     cdf = cdfplot(log10(pVal{i}));
%     cdf.LineStyle = linStyl;
%     cdf.LineWidth = 2;
%     cdf.Color = clr;
% 
% 
% end
% 
% pChance = 0;
% k = nTrial*Ps; %number of hits
% n = nTrial; %number of trials
% p = Ps; %probability of stimuli
% for j = k:n
%     Pj = nchoosek(n,j) * (p^j) * ((1-p)^(n-j)) ;
%     pChance = sum([pChance,Pj]);
% end
% xline(log10(pChance),'g-','LineWidth',2)
% xline(log10(0.05),'g:','LineWidth',2)
% xline(log10(0.01),'g--','LineWidth',2)
% % xline(Ps,'g--','LineWidth',2);
% xlabel('log10(pVal)')
% ylabel('proportion')
% title([anaMode ' Bayesian Decoder Performance'])
% legend({['v1bf; n=' num2str(length(hitRate{1}))],['v1af; n=' num2str(length(hitRate{2}))],['pssbf; n=' num2str(length(hitRate{3}))],['pssaf; n=' num2str(length(hitRate{4}))],'chance hitrate','pVal=0.05','pVal=0.01'})
% fName = fullfile(dsFold,['bayesDecoder_cdf_' anaMode 'pVal.fig']);
% saveas(gcf,fName)
% 











