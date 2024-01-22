%code written by Brandon Nanfito
%violin plot
function [v] = violin(dist,bw,lims,norm,plt)

% plt = 1;
% kernel = normpdf(-3:0.2:3);
% bw = 0.2;
n = length(dist);
s = std(dist);
rng = [min(dist) max(dist)];
for dp = 1:n

    if isnan(dist(dp))
        continue
    end

    curVal = dist(dp);
%     gX(dp,:) = linspace(rng(1)-(diff(rng)*1),rng(2)+(diff(rng)*1),10000);
    gX(dp,:) = linspace(lims(1),lims(2),lims(3));
    g(dp,:) = gauss(gX(dp,:),curVal,s*bw);

end
v.gX = gX;
v.g = g;
v.kdeY = gX(1,:);
v.kdeX = sum(g);
if norm == 1
    v.kdeX = v.kdeX/max(sum(g));
end


if plt==1
figure;hold on
plot(vertcat(zeros(size(dist')),ones(size(dist'))),repmat(dist',2,1),'k');
for dp = 1:n
    plot(g(dp,:)*-1,gX(1,:),'b')
end
plot(v.kdeX,v.kdeY,'r');
plot(v.kdeX*-1,v.kdeY,'r');
% ylim([0 1])
end


end
