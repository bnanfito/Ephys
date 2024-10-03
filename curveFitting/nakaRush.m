% naka rushton function


function [x,y,cF,resnorm,residuals,aic,bic] = nakaRush(r,c,mod,p)

if mod == 1
    
    x = linspace(1,100);

    rM = max(r);

    err = @(pars) (rM*((c.^pars(1))./((c.^(pars(3)*pars(1)))+(pars(2)^(pars(3)*pars(1))))))-r;
    [params,resnorm,residuals] = lsqnonlin(err,[1 50 1]);
    n = params(1);
    cF = params(2);
    s = params(3);

    %nakaruston fit
    y = rM* ( (x.^n) ./ ((x.^(s*n))+(cF^(s*n))) ) ;
    
    %Find the contrast for the nr fit that is closest to 50% response (the modified NakaRush gives)
    dif = abs(y-0.5);
    cF = find(dif == min(dif))/1000;
    
    if p == 1
        figure; hold on
        plot(c,r,'b','LineWidth',2)
        plot(x,y,'k','LineWidth',2)
        plot([0 cF],[rM/2 rM/2],'r')
        plot([cF cF],[0 rM/2],'r')
        xlim([0 100])
    end
    
elseif mod == 0
    
    x = linspace(1,100);

    rM = max(r);
    
    err = @(pars) (rM*((c.^pars(1))./((c.^pars(1))+(pars(2)^pars(1)))))-r;
    [params,resnorm,residuals] = lsqnonlin(err,[1 50]);
    n = params(1);
    cF = params(2);

    y = rM * ( (x.^n) ./ ((x.^n)+(cF^n)) ) ;
    
    if p == 1
        figure; hold on
        plot(c,r,'b','LineWidth',2)
        plot(x,y,'k','LineWidth',2)
        plot([0 cF],[rM/2 rM/2],'r')
        plot([cF cF],[0 rM/2],'r')
        xlim([0 100])
    end
    
end

n = length(c);
p = length(params);
rms = sqrt(resnorm/n);
aic = (n*log(rms))+(2*p);
bic = (n*log(rms))+(log(n)*p);

end