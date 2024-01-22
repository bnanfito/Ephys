function [ci] = confInt(x)
% compute the confidence interval for the population mean of x

% input 'x' is a vector containing samples of a variable. if its a 2D
% matrix, each column will be treated as a group.
x = squeeze(x);
nDim = max(size(size(x)));
if nDim>2

    disp('too many dimensions in data input')
    ci = nan;

elseif nDim == 2

    nGroups = size(x,2);

elseif nDim == 1
    
    nGroups = 1;

end

ci = nan(2,nGroups);

for g = 1:nGroups

    gDat = x(:,g);

    n = length(gDat);
    xbar = mean(gDat);
    se = std(gDat)/sqrt(n);
    nu = n - 1;
    
    conf = 0.95;
    alpha = 1 - conf;
    pLo = alpha/2;
    pUp = 1 - alpha/2;
    
    crit = tinv([pLo pUp], nu);
    
    ci(:,g) = xbar + crit*se;

end



end