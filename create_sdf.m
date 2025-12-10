%sdf

%compute sdf for a given vector of spike times (+trial IDs)

function [sdf] = create_sdf(dat)

    sf = 200;%sampling frequency (samples per second)
    bins = -1:1/sf:2;
    kernel = normpdf(-2:4/(sf*0.1):2,0,1);
    for u = 1:height(dat)
        nTrials = numel(dat.fr(u).trialNum);
        spkTimes = dat.spkTimes{u}(1,:);
        trialIds = dat.spkTimes{u}(2,:);
%         condIds = 
    
        for t = 1:nTrials
            sCounts = histcounts(spkTimes(trialIds==t),bins);
            unique(sCounts)
            sdf(t,:,u) = conv(sCounts,kernel,'same');
        end
    
        figure;hold on
%         plot(bins(2:end),sdf(:,:,u)')
        patch([bins(2:end) fliplr(bins(2:end))],[ mean(sdf(:,:,u))-sem(sdf(:,:,u)) fliplr(mean(sdf(:,:,u))+sem(sdf(:,:,u))) ],'k','FaceAlpha',0.2)
        plot(bins(2:end),mean(sdf(:,:,u)))
    end

end