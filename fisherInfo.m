%fi_bn
% calulate the fisher information metric for a tuning curve input

function [fi] = fisherInfo(r,c)

    %INPUT:
    % r = nRep x nCond matrix of responses
    % c = nCond vector of stimulus parameter values

    if size(c,1) == 2
        c1L = length(unique(c(1,:)));
        c2L = length(unique(c(2,:)));
        if c1L>c2L
            c = c(1,:);
        elseif c2L<c1L
            c = c(2,:);
        end
    end

    rMean = mean(r,'omitnan');
    rVar = std(r,'omitnan').^2;
    rSE = std(r,'omitnan')/sqrt(size(r,1));

    % First interpolate:
    cInterp = min(c):max(c);
    rMeanSmoothed = interp1(c,rMean,cInterp,'spline');
    rVarSmoothed  = interp1(c,rVar,cInterp,'spline');
    rVarSmoothed(rVarSmoothed<0) = 0;
    
    % get max slope (absolute value) and variance at nearest measured value
    absSlope = abs(gradient(rMeanSmoothed));
    FIslope = max(absSlope);
    FIvar = rVarSmoothed(absSlope==max(absSlope));
    
    % Fisher Info
    if FIslope == 0 || length(FIvar)>1
        fi = nan;
    else
        fi = FIslope^2 / FIvar;
    end
    % units are kinda arbitrary; if you did this for a bunch of neurons you
    % could just normalize the values to a range of 0-1 so it is comparable to
    % the other indices

%     figure;
%     subplot(2,2,1);hold on
%     errorbar(c,rMean,rSE)
% 
%     subplot(2,2,2);hold on
%     plot(cInterp,rMeanSmoothed)
% 
%     subplot(2,2,4);hold on
%     plot(absSlope)


end 