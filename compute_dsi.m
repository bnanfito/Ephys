%compute direction selectivity index (dsi) given vector of responses and
%directions

%inputs:
%   c = vector of direction (deg) corresponding to each column of r
%   r = matrix of responses size Reps x Conditions

function [dsi,rPref,cPref,dsi_rep,rP,cP,dsi_s,rP_s,cP_s,dsi_srep] = compute_dsi(c,r,plt)

        rMean = mean(r,1,'omitnan');
        rSem = std(r,[],1,'omitnan')./sqrt(size(r,1));
        rPref = max(rMean);
        if sum(rMean==rPref)>1
            rIn = rMean;
            pks = rIn == rPref;
            rConv = rIn([end 1:end 1]);
            rConv = conv(rConv,ones(1,3)*(1/3),'same');
            rConv = rConv(1+1:end-1);
            rConv(~pks) = 0;
            cPref = c(find(rConv==max(rConv),1,'first'));
            clear rIn pks rConv 
        else
            cPref = c(rMean==rPref);
        end
        cNull = mod(cPref+180,360);
        rNull = rMean(c==cNull);
        dsi = 1-(rNull/rPref);

        for rep = 1:size(r,1)
            rRep = r(rep,:);
            rP(rep) = max(rRep);
            if sum(rRep==rP(rep))
                rIn = rRep;
                pks = rIn == rP(rep);
                rConv = rIn([end 1:end 1]);
                rConv = conv(rConv,ones(1,3)*(1/3),'same');
                rConv = rConv(1+1:end-1);
                rConv(~pks) = 0;
                cP(rep) = c(find(rConv==max(rConv),1,'first'));
                clear rIn pks rConv 
            else
                cP(rep) = c(rRep==rP(rep));
            end
            cN(rep) = mod(cP(rep)+180,360);
            rN(rep) = rRep(c==cN(rep));
            dsi_rep(rep) = 1-(rN(rep)/rP(rep));


            wHan=hanning(3);
            wHan=wHan/sum(wHan);
            rRep_wrap = [rRep rRep rRep];
            rSRep_tmp = conv(rRep_wrap,wHan,'same');
            rSRep(rep,:) = rSRep_tmp((1:length(rRep))+length(rRep));
            rP_srep(rep) = max(rSRep(rep,:));
            if sum(rSRep(rep,:)==rP_srep(rep))>1
                rIn = rSRep(rep,:);
                pks = rIn == rP_srep(rep);
                rConv = rIn([end 1:end 1]);
                rConv = conv(rConv,ones(1,3)*(1/3),'same');
                rConv = rConv(1+1:end-1);
                rConv(~pks) = 0;
                cP_srep(rep) = c(find(rConv==max(rConv),1,'first'));
                clear rIn pks rConv 
            else
                cP_srep(rep) = c(rSRep(rep,:)==rP_srep(rep));
            end
            cN_srep(rep) = mod(cP_srep(rep)+180,360);
            rN_srep(rep) = rSRep(rep,c==cN_srep(rep));
            dsi_srep(rep) = 1-(rN_srep(rep)/rP_srep(rep));

        end

        wHan=hanning(3);
        wHan=wHan/sum(wHan);
        rWrap = [rMean rMean rMean];
        rSmooth = conv(rWrap,wHan,'same');
        rSmooth = rSmooth((1:length(rMean))+length(rMean));
        rP_s = max(rSmooth);
        if sum(rSmooth==rP_s)>1
            rIn = rSmooth;
            pks = rIn == rP_s;
            rConv = rIn([end 1:end 1]);
            rConv = conv(rConv,ones(1,3)*(1/3),'same');
            rConv = rConv(1+1:end-1);
            rConv(~pks) = 0;
            cP_s = c(find(rConv==max(rConv),1,'first'));
            clear rIn pks rConv 
        else
            cP_s = c(rSmooth==rP_s);
        end
        cN_s = mod(cP_s+180,360);
        rN_s = rSmooth(c==cN_s);
        dsi_s = 1-(rN_s/rP_s);
        

        if plt==1
            figure; 
            subplot(1,2,1);hold on
            yline(0,'k')

            plot(c,r','c')
            plot(cP,rP,'co','MarkerFaceColor','c')
            plot(cN,rN,'co','MarkerFaceColor','none')

            errorbar(c,rMean,rSem,'b','LineWidth',2)
            plot(cPref,rPref,'bo','MarkerFaceColor','b','MarkerSize',10)
            plot(cNull,rNull,'bo','MarkerFaceColor','none','MarkerSize',10,'LineWidth',2)
            yline(rPref,'b')
            yline(rNull,'b--')
%             plot(repmat(mean([cPref,cNull]),1,2),[rPref,rNull],'b','LineWidth',2)
%             plot(repmat(mean([cPref,cNull]),1,2),[0,rNull],'b--','LineWidth',2)

            plot(c,rSRep','m')
            plot(cP_srep,rP_srep,'mo','MarkerFaceColor','m')
            plot(cN_srep,rN_srep,'mo','MarkerFaceColor','none')

            plot(c,rSmooth,'r','LineWidth',2)
            plot(cP_s,rP_s,'ro','MarkerFaceColor','r','MarkerSize',10)
            plot(cN_s,rN_s,'ro','MarkerFaceColor','none','MarkerSize',10,'LineWidth',2)
            yline(rP_s,'r')
            yline(rN_s,'r--')

            subplot(1,2,2);hold on
            plot(1:size(r,1),dsi_rep,'co','MarkerFaceColor','c')
            yline(mean(dsi_rep),'c','LineWidth',2)
            yline(dsi,'b','LineWidth',2)
            plot(1:size(r,1),dsi_srep,'mo','MarkerFaceColor','m')
            yline(mean(dsi_srep),'m','LineWidth',2)
            yline(dsi_s,'r','LineWidth',2)
            ylim([0 1])
        end

end