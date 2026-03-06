function [vm] = fit_vonMises(r,t)
t = deg2rad(t);
plr=0;

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

thetaFine = linspace(0, 2*pi, 360);

Rfit = doubleVonMises(params_fit, thetaFine);

vm.rho = Rfit;
vm.theta = rad2deg(thetaFine);
vm.params = params_fit;
vm.Rp = vm.params(1) + vm.params(2)*exp(vm.params(4)) + vm.params(3)*exp(-vm.params(5));
vm.Rn = vm.params(1) + vm.params(3)*exp(vm.params(5)) + vm.params(2)*exp(-vm.params(4));
vm.Tp = rad2deg(params_fit(end));

% figure
% if plr == 1
% polarplot(t, r, 'ko','MarkerFaceColor','k')
% hold on
% polarplot(deg2rad(vm.theta), vm.rho, 'r','LineWidth',2)
% else
% plot(t, r, 'ko','MarkerFaceColor','k')
% hold on
% plot(deg2rad(vm.theta), vm.rho, 'r','LineWidth',2)
% end
% title('Double von Mises Fit')

end