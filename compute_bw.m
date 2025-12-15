%compute_bw

function [BW,BWSmooth] = compute_bw(c,r)

    [prefResp,prefCond]=max(r);

    if length(find(r==prefResp))==1
        %only one max - can use it directly
        prefDir=c(prefCond);
    else
        %use convolution to pick the maximum that is surrounded by
        %higher values; using square convolution   
        tc=[r(end) r r(1)];  %wrap around tuning curve first (need 1 on either side)
        idx=find(r==prefResp)+1; %indices for peaks
        rateConv=conv(tc,[1 1 1],'same'); %convolution
        mxConv=rateConv(idx); %these are the convolution values for the peaks, pick the bigger one    
        mx=idx(mxConv==max(mxConv)); 
        mx=mx(1);%if there are multiple, we'll just take the first
        prefDir=c(mx-1);
    end
    
    %% preferred response
    prefResp=r(c==prefDir);
    
    %% null direction and response
    nullDir=mod(prefDir+180,360);
    nullResp=r(c==nullDir);
    
    %% OSI
    DSI=(prefResp-nullResp)/prefResp;
    
    %% vector length - better with rectified response
    avgRateR=r;
    avgRateR(avgRateR<0)=0;
    dirVec=sum(avgRateR.*exp(1i*c*pi/180)); 
    Ldir=abs(dirVec)/sum(avgRateR);
    Adir=mod(angle(dirVec)*180/pi/2,180);
    
    %% bandwidth - as in ringach 2002, with and without smoothing
    wHan=hanning(3);
    wHan=wHan/sum(wHan);
    % avgRateS=conv(avgRate,wHan,'same');
    
    %wrap around to avoid edge arfecats (need to do that before interpolation
    tcWrap=[r r r];
    % tcWrapS=[avgRateS avgRateS avgRateS];
    tcWrapS=conv(tcWrap,wHan,'same');
    avgRateS=tcWrapS((1:length(c))+length(c));
    
    %recompute the preferred direction (and associated metrics) for smoothed tc
    %because smoothing can change the relative amplitude of peaks (ie change pref)
    [prefRespS,prefCondS] = max(avgRateS);
    if length(find(avgRateS==prefRespS))==1
        prefDirS=c(prefCondS);
    else
        tc=[avgRateS(end) avgRateS avgRateS(1)];
        idx=find(avgRateS==prefRespS)+1;
        rateConv=conv(tc,[1 1 1],'same');
        mxConv=rateConv(idx);
        mx=idx(mxConv==max(mxConv));
        mx=mx(1);
        prefDirS=c(mx-1);
    end
    prefRespS=avgRateS(c==prefDirS);
    nullDirS=mod(prefDirS+180,360);
    nullRespS=avgRateS(c==nullDirS);
    DSI_S=(prefRespS-nullRespS)/prefRespS;
    
    %criterion response: peak resp/sqrt(2)
    critResp=prefResp/sqrt(2);
    critRespSmooth=prefRespS/sqrt(2);
    
    %interpolate tuning function and orientation vector
    nInterp=3000; %3*1000 because of wrap around
    domInter=linspace(c(1)-360,c(end)+360,nInterp); %equivalent to 3x range
    
    tcIP=griddedInterpolant(tcWrap); %generates object to be evaluated
    tcInter=tcIP(linspace(1,length(tcWrap),nInterp));
    
    tcIPS=griddedInterpolant(tcWrapS); %smooth version
    tcInterS=tcIPS(linspace(1,length(tcWrapS),nInterp));
    
    %find crossings with criterion response
    [~,~,cIdx]=zerocrossrate(tcInter,'Level',critResp);
    cIdx(1)=0;
    if isempty(cIdx) || sum(cIdx)<4 %less than 4: only 1 data point of the tc falls below the level
        cross1=nan;
        cross2=nan;
        BW=NaN;
    else
        %find the ones closest to the peak
        crossOri=domInter(cIdx);
        cross1=max(crossOri(crossOri<prefDir));
        cross2=min(crossOri(crossOri>prefDir));
        BW=(cross2-cross1)/2;
    end
    
    [~,~,cIdxS]=zerocrossrate(tcInterS,'Level',critRespSmooth);
    cIdxS(1)=0;
    if isempty(cIdxS) || sum(cIdxS)<4
        cross1S=nan;
        cross2S=nan;
        BWSmooth=NaN;
    else
        %find the ones closest to the peak
        crossOri=domInter(cIdxS);
        cross1S=max(crossOri(crossOri<prefDirS));
        cross2S=min(crossOri(crossOri>prefDirS));
        BWSmooth=(cross2S-cross1S)/2;
    end


end

