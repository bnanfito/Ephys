
function [g] = dirGauss(r,dir,p)

X = [-360:720];
r_rep = repmat(r,3,1);
dir_rep = [dir-360;dir;dir+360];

%% shift to minimum

Rmin = min(r);
dirMin = dir(r == Rmin);
if length(dirMin)>1
%     diffs = abs(dirMin-(dir(r==max(r))))
%     dirMin = dirMin(diffs == min(diffs));
    dirMin = min(dirMin);
end
if length(dirMin)>1
    dirMin = dirMin(1);
end
dirShift = [dir(dir<dirMin)+360 dir(dir>=dirMin)];
if mod(max(dirShift),360) == dir(r==max(r)) & ~isempty(dir(dir<dirMin))
    dirShift(find(dirShift == max(dirShift))+1)=dirShift(find(dirShift == max(dirShift))+1)+360;
end

if p
    figure;
    subplot(2,2,1);hold on
    plot(dir,r,'ko')
    plot(dirShift,r,'r*')
    plot(dir_rep,r_rep,'k-')
    plot([0 0],[0 max(r)],'r--')
    plot([360 360],[0 max(r)],'r--')
end

%% initial guess at gaussian

Rpref = max(r);
pref = unique(dirShift(r == Rpref));
if length(pref)>1
    for i = 1:length(pref)
        id = find(dirShift == pref(i));
        if length(id)>1 && (ismember(1,id) && ismember(length(dirShift),id))
            energies(i) = mean([r(end-1) r(1) r(2)]);
        else
            energies(i) = mean(r(id-1:id+1));
        end
    end
    pref = pref(energies == max(energies));
    if length(pref)>1
        pref = pref(1);
    end
end
sigma1 = 30;
nullA = pref+180;
nullB = pref-180;
if ismember(nullA,dirShift)
null = nullA;
elseif ismember(nullB,dirShift)
null = nullB;
end
Rnull = r(dirShift == null);
sigma2 = 30;

gauss = @(x) Rpref*exp(-((x-pref).^2/(2*sigma1^2))) + Rnull*exp(-((x-null).^2/(2*sigma2^2)));


if p
    subplot(2,2,2);hold on
    plot(dir,r,'ko')
    plot(dirShift,r,'r*')
    plot(pref,Rpref,'ro')
    plot(null,Rnull,'ro')
    plot(X,gauss(X),'g--')
    plot([0 0],[0 max(r)],'r--')
    plot([360 360],[0 max(r)],'r--')
end

%% fitting gaussian

err = @(pars) (pars(1)*exp(-((dirShift-pars(2)).^2/(2*pars(3)^2))) + pars(4)*exp(-((dirShift-pars(5)).^2/(2*pars(6)^2))))-r;
% err = @(pars) (pars(1)*exp(-(mod(dirShift-pars(2),360).^2/(2*pars(3)^2))) + pars(4)*exp(-(mod(dirShift-pars(5),360).^2/(2*pars(6)^2))))-r;

lb = [0     pref-45     10      0       null-45     10];
ub = [Inf   pref+45     100     Inf     null+45     100];
xo = [Rpref , pref , sigma1 , Rnull , null, sigma2];
[params,resnorm,residuals] = lsqnonlin(err,xo,lb,ub);

Rpref = params(1);%max(r);
pref = params(2);%dir(r == max(r));
sigma1 = params(3);%10;
Rnull = params(4);%r(dir == null);
null = params(5);%mod(pref+180,360);
sigma2 = params(6);%10;

g.auss = @(x) Rpref*exp(-((x-pref).^2/(2*sigma1^2))) + Rnull*exp(-((x-null).^2/(2*sigma2^2)))...
    + Rnull*exp(-((x-(null-360)).^2/(2*sigma2^2))) + Rpref*exp(-((x-(pref-360)).^2/(2*sigma1^2)));
g.Rpref = Rpref;
g.pref = pref;
g.sigma1 = sigma1;
g.Rnull = Rnull;
g.null = null;
g.sigma2 = sigma2;
g.residuals = g.auss(dir)-r;


% g.auss = @(x) Rpref*exp(-(mod(x-pref,360).^2/(2*sigma1^2))) + Rnull*exp(-(mod(x-null,360).^2/(2*sigma2^2)));
% g.auss = @(x) Rpref*exp(-((mod(x,360)-mod(pref,360)).^2/(2*sigma1^2))) + Rnull*exp(-((mod(x,360)-mod(null,360)).^2/(2*sigma2^2)));


if p
    subplot(2,2,3);hold on
    plot(dir,r,'ko')
    plot(dirShift,r,'r*')
    plot(X,g.auss(X),'g-')
    plot([0 0],[0 max(r)],'r--')
    plot([360 360],[0 max(r)],'r--')
end













% minR = min(r);
% minDir = dir(r == minR);
% if length(minDir)>1
%     minDir = minDir(1);
% end
% shift = -(find(dir == minDir)-1);
% 
% r = circshift(r, shift)';
% dir = circshift(dir, shift)';
% length(dir)+(shift+1):length(dir);
% dir(length(dir)+(shift+1):length(dir)) = dir(length(dir)+(shift+1):length(dir))+360;
% 
% if p
%     figure;hold on
%     plot(mod(dir,360),r,'ko','MarkerSize',10)
% end
% 
% X = linspace(0,359,360);
% 
% Rpref = max(r);
% pref = dir(r == max(r));
% if ismember(pref+180,dir)
%     null = pref+180;
% elseif ismember(pref-180,dir)
%     null = pref-180;
% end
% Rnull = r(mod(dir,360) == mod(null,360));
% sigma1 = 10;
% sigma2 = 10;
% 
% if p
% initial = Rpref*exp(-((X-mod(pref,360)).^2/(2*sigma1^2))) + Rnull*exp(-((X-mod(null,360)).^2/(2*sigma2^2)));
% plot(mod(X,360),initial,'g','LineWidth',3)
% end
% 
% err = @(pars) (pars(1)*exp(-((dir-pars(2)).^2/(2*pars(5)^2))) + pars(4)*exp(-((dir-pars(3)).^2/(2*pars(6)^2))))-r;
% 
% lb = [0 pref-45 null-45 0 0 0];
% ub = [Inf pref+45 null+45 Inf 45 45];
% xo = [Rpref , pref , null , Rnull , sigma1, sigma2];
% [params,resnorm,residuals] = lsqnonlin(err,xo,lb,ub);
% 
% Rpref = params(1);%max(r);
% pref = mod(params(2),360);%dir(r == max(r));
% null = mod(params(3),360);%mod(pref+180,360);
% Rnull = params(4);%r(dir == null);
% sigma1 = params(5);%10;
% sigma2 = params(6);%10;
% 
% g.auss = @(x) Rpref*exp(-((x-pref).^2/(2*sigma1^2))) + Rnull*exp(-((x-null).^2/(2*sigma2^2)));
% g.auss1 = @(x) Rpref*exp(-((x-pref).^2/(2*sigma1^2)));
% g.auss2 = @(x) Rnull*exp(-((x-null).^2/(2*sigma2^2)));
% g.Rpref = Rpref;
% g.pref = pref;
% g.Rnull = Rnull;
% g.null = null;
% g.sigma1 = sigma1;
% g.sigma2 = sigma2;
% 
% if p
% plot(X,g.auss(X) + g.auss(X-360) + g.auss(X+360),'k','LineWidth',2)
% plot(X,g.auss1(X) + g.auss1(X-360) + g.auss1(X+360),'b--')
% plot(X,g.auss2(X) + g.auss2(X-360) + g.auss2(X+360),'r--')
% 
% axis tight
% 
% % labels = mod(xticks,360);
% % xticklabels(labels)
% end





%     x = linspace(0,359,360);
% 
%     err = @(pars) (pars(1) + pars(2)*exp(-((dir-pars(3)).^2/(2*pars(5)^2))) + pars(4)*exp(-((dir-(pars(3)+180)).^2/(2*pars(6)^2))))-r;
% 
%     lb = [-Inf 0 0 0 -Inf -Inf];
%     ub = [Inf max(r) 359 max(r) Inf Inf];
%     xo = [min(r) , max(r)-min(r) , dir(r==max(r)) , r(dir==mod(dir(r==max(r))+180,360))-min(r) , 10, 10];
%     [params,resnorm,residuals] = lsqnonlin(err,xo,lb,ub);
% 
%     Roff = params(1);%min(r);
%     Rpref = params(2);%max(r)-Roff;
%     pref = params(3);%dir(r == max(r));
%     null = pref + 180;
%     Rnull = params(4);%r(dir == null)-Roff;
%     sigma1 = params(5);%10;
%     sigma2 = params(6);%10
%     
%     g.auss = Roff + Rpref*exp(-((x-pref).^2/(2*sigma1^2))) + Rnull*exp(-((x-null).^2/(2*sigma2^2)));
%     g.Roff = Roff;
%     g.Rpref = Rpref;
%     g.pref = pref;
%     g.Rnull = Rnull;
%     g.null = null;
%     g.sigma1 = sigma1;
%     g.sigma2 = sigma2;

    
    
    
    
    
    
    
%     figure;hold on
%     plot(dir, r, 'ko')
%     plot(x,g,'r')

%     err = @(pars) (pars(1)*exp(((-0.5*(dir-pars(2)).^2)/pars(3)^2)+pars(4)))-r;
%     [params,resnorm,residuals] = lsqnonlin(err,[2,2,2,2]);
%     
%     a = params(1);
%     b = params(2);
%     c = params(3);
%     d = params(4);
% 
%     g = a*exp(((-0.5*(x-b).^2)/c^2)+d);
%     
%     figure;hold on
%     plot(dir,r,'ko')
%     plot(x,g)

end