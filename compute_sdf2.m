function [sdf,time] = compute_sdf2(spkTimes,dt,time)

w = dt/10;
g = zeros(1,length(time));
for s = 1:length(spkTimes)
    g(s,:) = normpdf(time,spkTimes(s),w)*dt;
end
sdf = sum(g,1)/dt;

end