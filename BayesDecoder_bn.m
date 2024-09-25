%Bayesian Decoder
%written by Brandon Nanfito

clear all
close all

anaMode = 'SU';
load(['/Volumes/Lab drive/Brandon/data/dataSets/training/Train_V1Cool/' anaMode '/Train_V1Cool_' anaMode 'dataSet.mat'])

figure;hold on
for i = 1:4

    if i==1
        tbl = data.v1bf;
        clr = 'b';
        linStyl = '--';
    elseif i==2
        tbl = data.v1af;
        clr = 'b';
        linStyl = '-';
    elseif i==3
        tbl = data.pssbf;
        clr = 'r';
        linStyl = '--';
    elseif i==4
        tbl = data.pssaf;
        clr = 'r';
        linStyl = '-';
    end
    
    nU = height(tbl);
    for u = 1:nU
    
        if ~tbl.goodUnit(u)
            continue
        end
    
        R = tbl.response{u};
        rVar{i}(u) = mean(var(R,'omitnan'));
        nStim = size(R,2);
        nRep = size(R,1);
        nTrial = nStim*nRep;
        R = R(:);
        if strcmp(anaMode,'MU')
            skipTrial = isnan(R);
        end
        S = repmat(1:nStim,nRep,1);
        S = S(:);
        for t = 1:nTrial

            if strcmp(anaMode,'MU') && ismember(t,find(skipTrial))
                continue
            end
    
            Pr = sum(R==R(t))/nTrial;
            for s = 1:nStim
                Ps = 1/nStim;
                Pr_s = sum(R(S==s)==R(t))/nRep;
    
                Ps_r{i}(s,t,u) = (Pr_s*Ps)/(Pr);
            end
    
            guess = find(Ps_r{i}(:,t,u)==max(Ps_r{i}(:,t,u)));
            if length(guess)>1
                guess = guess(randi(length(guess),1));
            end
            sP{i}(t,u) = guess;
    
        end
    
        if strcmp(anaMode,'MU')
            perf{i}(u) = sum(S(~skipTrial)==sP{i}(~skipTrial,u))/nTrial;
        else
            perf{i}(u) = sum(S==sP{i}(:,u))/nTrial;
        end
    
    %     figure
    %         subplot(2,2,1)
    %         histogram(R(:))
    % 
    %         for s = 1:nStim
    %             subplot(nStim,2,2*s);hold on
    %             histogram(R(S==s))
    %         end
    % 
    %         title(['unit# ' num2str(u)])
    
    end
    
    perf{i} = perf{i}(tbl.goodUnit);
    rVar{i} = rVar{i}(tbl.goodUnit);
    cdf = cdfplot(perf{i});
    cdf.LineStyle = linStyl;
    cdf.Color = clr;


end

xline(Ps,'g--')


