%compute direction selectivity index (dsi) given vector of responses and
%directions

function [dsi,rPref,cPref] = compute_dsi(c,r)

        rPref = max(r);
        if sum(r==rPref)>1
            rIn = r;
            pks = rIn == rPref;
            rConv = rIn([end 1:end 1]);
            rConv = conv(rConv,ones(1,3)*(1/3),'same');
            rConv = rConv(1+1:end-1);
            rConv(~pks) = 0;
            cPref = c(find(rConv==max(rConv),1,'first'));
            clear rIn pks rConv 
        else
            cPref = c(r==rPref);
        end
        cNull = mod(cPref+180,360);
        rNull = r(c==cNull);
        dsi = 1-(rNull/rPref);

end