
function [g] = gauss(x,m,s)
% a = mean
% sigma = std
% x = evaluation

% g = exp(-1*( ( 2*(x-m) )/s ).^2);

g = 1/sqrt(2*pi*(s^2)) * exp(-1*( ((x-m).^2)/(2*(s^2)) ));

end