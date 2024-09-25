%decoder performance (p-value for binary hit rates)
clear all
close all

K = 16; % number of possible stimuli
n = 80;  % number of trials
Ps = 1/K; % probability of a given stimuli
for k = 1:n % number of hits

pVal(k) = 0;
for j = k:n
    Pj = nchoosek(n,j) * (Ps^j) * ((1-Ps)^(n-j));
    pVal(k) = sum([pVal(k),Pj]);
end

end

figure; plot((1:n)/n,log10(pVal))
xlabel('hit rate')
ylabel('log10(pVal)')