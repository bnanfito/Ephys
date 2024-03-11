
function [g] = gauss(x,a,b)
% a = mean
% b = std
% x = evaluation

g = exp(-1*( ( 2*(x-a) )/b ).^2);

end