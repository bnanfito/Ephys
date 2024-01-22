
function [mv] = meanvec(theta,rho)

mv.com = sum(rho(:).*exp(sqrt(-1)*mod(deg2rad(theta(:)),2*pi))); %resultant (mean) vector; the weighted sum of cos & sin of angles in radians
mv.mag = abs(mv.com); %resultant vector magnitude 
mv.ang = rad2deg(mod(angle(mv.com),2*pi)); %resultant vector angle in deg
mv.cv = 1-((mv.mag)./sum(rho(:))); %circular variance; 1 minus the normalized magnitude of the resultant vector (the magnitude of the weighted sum of angles normalized by the sum of responses)


end