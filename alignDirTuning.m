
function [tOut,rOut,i] = alignDirTuning(tIn,rIn)


    rPref = max(rIn);
    if length(find(rIn == rPref))>1

        pks = rIn == rPref;
        rConv = rIn([end 1:end 1]);
        rConv = conv(rConv,ones(1,3)*(1/3),'same');
        rConv = rConv(1+1:end-1);
        rConv(~pks) = 0;
        cPref = find(rConv==max(rConv),1,'first');

    else

        cPref = find(rIn == rPref);

    end

    tCent = tIn - tIn(cPref);

    tCent(tCent>180) = tCent(tCent>180)-360;
    tCent(tCent<-180) = tCent(tCent<-180)+360;

    [tOut,i] = sort(tCent);
    rOut = rIn(i);
    if abs(tOut(1)) < abs(tOut(end))
        tOut = [tOut(end)*-1 tOut];
        i = [i(end) i];
        rOut = [rOut(end) rOut];
    elseif abs(tOut(1)) > abs(tOut(end))
        tOut = [tOut tOut(1)*-1];
        i = [i i(1)];
        rOut = [rOut rOut(1)];
    end



end