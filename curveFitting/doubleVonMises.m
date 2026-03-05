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