%information theory

clear all;close all
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
nRep = 10;
% s = 0:22.5:337.5;
% rMean = normpdf(s,180,40)*500; 
% rStd = repmat(1,1,length(rMean));
s = 0:10;
rMean = [20 20 20 30 40 50 60 70 80 80 80];
rStd = repmat(10,1,11);
r_sim = rMean+(randn(nRep,length(s)).*rStd); 
rMean_sim = mean(r_sim,'omitnan');
rVar = rStd.^2;
rStd_sim = std(r_sim,'omitnan');
rVar_sim = rStd_sim.^2;

% R_axMin = min(r_sim,[],'all');R_axMax = max(r_sim,[],'all');
R_axMin = 0;R_axMax = 100;
R_ax = R_axMin:1:R_axMax; 
R_lbl = R_axMin:2:R_axMax; 
R_tic = find(ismember(R_ax,R_lbl));

pR_S = normpdf(repmat(R_ax,length(s),1)',rMean,rStd);
pR = sum(pR_S,2)/size(pR_S,2); 
pRsim_S = normpdf(repmat(R_ax,length(s),1)',rMean_sim,rStd_sim);
pRsim = sum(pRsim_S,2)/size(pRsim_S,2);

figure(1)
subplot(1,2,1)
hold on
plot(s,rMean,'-o','Color',colors{1}); 
plot(repmat(s,nRep,1)',r_sim','.','Color',colors{2});
plot(s,rMean_sim,'-o','Color',colors{2})
xticks(s);
xlabel('ori')
ylabel('response (Hz)')
axis tight

figure(2)
subplot(2,2,1)
imagesc(pR_S)
set(gca,'YDir','normal')
title('P(r|s)')
yticks(R_tic)
yticklabels(R_lbl)

subplot(2,2,2)
plot(pR,R_ax)

subplot(2,2,3)
imagesc(pRsim_S)
set(gca,'YDir','normal')
title('P(r)')
yticks(R_tic)
yticklabels(R_lbl)

subplot(2,2,4)
plot(pRsim,R_ax)

for x = 1:length(s)
    Ir_s(x) = sum(pR_S(:,x).*log2(pR_S(:,x)./pR));
    Ir_s_sim(x) = sum(pRsim_S(:,x).*log2(pRsim_S(:,x)./pRsim));
end

figure(1)
subplot(1,2,2)
hold on
plot(s,Ir_s,'-o')
plot(s,Ir_s_sim,'-o')
xticks(s)
xlabel('ori')
ylabel('info (bits)')
axis tight