
function [g] = gauss(x,a,b)

g = exp(-1*( ( 2*(x-a) )/b ).^2);

end