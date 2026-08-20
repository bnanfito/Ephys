function [sdf,time] = compute_sdf2(spkTimes,dt,tLims)

time = tLims(1):dt:tLims(2);
w = 0.002;
g = nan(1,length(time));
for s = 1:length(spkTimes)
    g(s,:) = normpdf(time,spkTimes(s),w)*dt;
end
sdf = sum(g,1)/dt;

end