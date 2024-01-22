%convertAL

function sumStats = convertAL(animal,unit,expt,probe,dataFold)

alignBit = 1;

exptName = [animal '_u' unit '_' expt];
physDir = fullfile(dataFold,'AL',animal,exptName);
if strcmp(animal,'FEAQ7')
    load(fullfile(physDir,[exptName '_data.mat']))
else
    load(fullfile(physDir,[exptName '_p' num2str(probe) '_data.mat']))
end
data = Data{1,1}; clear Data
nUnits = size(UnitType{1,1},1)

c = CondInfo;
r = data.RespMean';
for u = 1:nUnits

    rPref(u,1) = max(r(:,u));
    pks = r(:,u)==rPref(u,1);
    if sum(pks)>1
        rConv = r(:,u);
        rConv = rConv([end 1:end 1]);
        rConv = conv(rConv,ones(1,3)*(1/3),'same');
        rConv = rConv(1+1:end-1);
        rConv(~pks) = 0;
        cPref(u,1) = c(find(rConv==max(rConv),1,'first'));
    else
        cPref(u,1) = c(r(:,u)==rPref(u,1));
    end
    cNull(u,1) = mod(cPref(u,1)+180,360);
    rNull(u,1) = r(c==cNull(u,1),u);

    dsi(u,1) = abs(rPref(u,1)-rNull(u,1))/rPref(u,1);
    mv = meanvec(c,r(:,u));
    dcv(u,1) = mv.cv;

    if alignBit == 1
        [tOut,rOut,i] = alignDirTuning(c',r(:,u)');
        tuningX{u,1} = tOut;
        tuningY{u,1} = rOut;
    else
        tuningX{u,1} = c;
        tuningY{u,1} = r(:,u);
    end

    isAct = rPref(u,1) > 2;
    isVis = anova1([data.AllBResp(:,u) squeeze(data.AllResp(u,:,:))],[],"off")<0.05;
    goodUnit(u,1) = isAct&isVis;


end

varNames = {'goodUnit','dsi','dcv','rPref','cPref','tuningX','tuningY'};
sumStats = table(goodUnit,dsi,dcv,rPref,cPref,tuningX,tuningY,'VariableNames',varNames);

end