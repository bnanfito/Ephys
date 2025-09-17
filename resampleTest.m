
function p = resampleTest(d1,d2)
% ag = 4;
% d1 = D.latency(idx3&D.AG==ag); d1 = d1(~isnan(d1));
% d2 = D.latency(idx4&D.AG==ag); d2 = d2(~isnan(d2));

if length(d1)<length(d2)
    x = d1;
    y = d2;
elseif length(d2)<length(d1)
    x = d2;
    y = d1;
end

for i = 1:10000
    Y(i,:) = datasample(y,length(x),'Replace',true);
    d(i) = median(y)-median(Y(i,:));
end

figure; tiledlayout(2,2)

nexttile; hold on
cX = cdfplot(x);
cX.Color = 'r';
xline(median(x),'r--')
cY = cdfplot(y);
cY.Color = 'b';
xline(median(y),'b--')
legend([cX,cY],{['n=' num2str(length(x))],['n=' num2str(length(y))]})

nexttile;hold on
cX = cdfplot(x);
cX.Color = 'r';
cX.LineWidth = 2;
xline(median(x),'r--','LineWidth',2)
cY = cdfplot(y);
cY.Color = 'b';
cY.LineWidth = 2;
xline(median(y),'b--','LineWidth',2)
for i = 1:100
    cY = cdfplot(Y(i,:));
    cY.Color = 'b';
    xline(median(Y(i,:)),'b--')
end

nexttile; hold on
histogram(d)
xline(median(x)-median(y))

p = [];

end