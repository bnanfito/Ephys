function [thetaFine,Rfit] = fit_vonMises(r,t)

% Initial parameter guesses
R0_init     = min(r);
A1_init     = max(r) - min(r);
A2_init     = A1_init / 2;
k1_init  = 2;
k2_init  = 2;
theta0_init = t(r == max(r));
theta0_init = theta0_init(1);

params0 = [R0_init, A1_init, A2_init, k1_init, k2_init, theta0_init];

% Bounds
lb = [0, 0, 0, 0, 0, -pi];
ub = [Inf, Inf, Inf, 100, 100, pi];

% Fit
options = optimoptions('lsqcurvefit','Display','off');
params_fit = lsqcurvefit(@doubleVonMises, params0, t, r, lb, ub, options);

thetaFine = linspace(0, 2*pi, 500);

Rfit = doubleVonMises(params_fit, thetaFine);

% figure
% polarplot(t, r, 'ko','MarkerFaceColor','k')
% hold on
% polarplot(thetaFine, Rfit, 'r','LineWidth',2)
% title('Double von Mises Fit')

end