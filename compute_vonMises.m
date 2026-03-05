function [thetaFine,Rfit] = compute_vonMises(r,t)

% theta = deg2rad(0:359);
% 
% R0 = 10;
% A1 = 5;
% A2 = 5;
% theta0 = deg2rad(180);
% kappa = 1;
% 
% R = R0 + (A1*exp(kappa * cos(theta-theta0))) + (A2*exp(kappa * cos(theta-theta0-pi)));
% figure;
% polarplot(theta,R)

function R = doubleVonMises(params, theta)

% params = [R0, A1, A2, theta0, kappa]

R0     = params(1);
A1     = params(2);
A2     = params(3);
k1     = params(4);
k2     = params(5);
theta0 = params(6);


% Ensure theta is column
theta = theta(:);
vm1 = A1 * exp(k1 * cos(theta - theta0));
vm2 = A2 * exp(k2 * cos(theta - theta0 - pi));
R = R0 + vm1 + vm2;

end

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
params_fit = lsqcurvefit(@doubleVonMises, ...
                         params0, ...
                         t, ...
                         r, ...
                         lb, ub, ...
                         options);

thetaFine = linspace(0, 2*pi, 500);

Rfit = doubleVonMises(params_fit, thetaFine);

% figure
% polarplot(t, r, 'ko','MarkerFaceColor','k')
% hold on
% polarplot(thetaFine, Rfit, 'r','LineWidth',2)
% title('Double von Mises Fit')

end