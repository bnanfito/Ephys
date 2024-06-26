%fitSig

function [params,residuals] = fitSig(inX,inY)

figure;hold on
plot(inX,inY,'-o')

err = @(pars) ( pars(1)./( 1+exp( -pars(2) * (inX-pars(3)) ) ) )-inY;

lb = [-inf 0.5 min(inX)];
ub = [inf inf max(inX)];
xo = [1 1 median(inX)];
[params,resnorm,residuals] = lsqnonlin(err,xo,lb,ub);

xSpace = min(inX):max(inX);

a = params(1);
b = params(2);
c = params(3);
sig = @(x) a./( 1+exp( -b * (x-c) ) );
plot(xSpace,sig(xSpace))


end







% clear all
% close all
% 
% load('/Users/brandonnanfito/Documents/NielsenLab/data/forAugusto/cntrl/DPI pre1 stim1/cntrl_dpi_data.mat')
% 
% for hr = 1:8
% blocks = [(hr*2)-1,hr*2];
% v1tbl{hr} = vertcat(v1{blocks(1),end},v1{blocks(2),end});
% v1tbl{hr} = v1tbl{hr}(v1tbl{hr}.goodUnit==1,:);
% v1dpiMean(hr) = mean(v1tbl{hr}.dpi);
% v1dpiMed(hr) = median(v1tbl{hr}.dpi);
% v1dpiSEM(hr) = std(v1tbl{hr}.dpi)/sqrt(length(v1tbl{hr}.dpi));
% if hr>1
% v1hrDiff(hr-1,:) = v1dpiMean()
% end
% psstbl = vertcat(pss{blocks(1),end},pss{blocks(2),end});
% psstbl = psstbl(psstbl.goodUnit==1,:);
% pssdpiMean(hr) = mean(psstbl.dpi);
% pssdpiMed(hr) = median(psstbl.dpi);
% pssdpiSEM(hr) = std(psstbl.dpi)/sqrt(length(psstbl.dpi));
% end
% 
% figure; 
% subplot(2,1,1);hold on
% y = v1dpiMean;
% sem = v1dpiSEM;
% x = 1:length(y);
% plot(y)
% plot(repmat(x,2,1),y+([1;-1]*sem))
% y = pssdpiMean;
% sem = pssdpiSEM;
% x = 1:length(y);
% plot(y)
% plot(repmat(x,2,1),y+([1;-1]*sem))
%  
% subplot(2,1,2);hold on
% yline(0,'k--')
% y = v1dpiMean(2:end)-v1dpiMean(1:end-1);
% plot(y)
% y = pssdpiMean(2:end)-pssdpiMean(1:end-1);
% plot(y)