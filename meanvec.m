
function [mv] = meanvec(theta,rho)

% INPUT
% 1. theta: vector containing the stimulus condition (direction of motion
% in deg) of each response in rho
% 2. rho: vector of responses

rho(rho<0) = 0; % rectify rho (ignore negative BCFR in vector space)

%% Dir vector space
% resultant (mean) vector; the weighted sum of cos & sin of angles in radians 
    mv.resultant_dir = sum(rho(:).*exp(sqrt(-1)*mod(deg2rad(theta(:)),2*pi))); 
% resultant vector angle in deg
    mv.angDir = rad2deg(mod( angle(mv.resultant_dir) ,2*pi)); 
% resultant vector magnitude
    mv.magDir = abs( mv.resultant_dir ); 
% normalized vector length in dir space (Ldir)
    mv.ldr = abs( mv.resultant_dir / sum(rho(:)) ); 
% directional circular variance
    mv.dcv = 1-mv.ldr;

%% Ori vector space
% resultant (mean) vector; the weighted sum of cos & sin of angles in radians 
    mv.resultant_ori = sum( rho(:) .* exp(sqrt(-1)*2*mod(deg2rad(theta(:)),pi)) ); 
% resultant vector angle in deg
    mv.angOri = rad2deg(mod( angle(mv.resultant_ori)/2 ,pi));
% resultant vector magnitude
    mv.magOri = abs(mv.resultant_ori); 
% normalized vector length
    mv.lor = abs( mv.resultant_ori / sum(rho(:)) );
% orientation circular variance
    mv.ocv = 1-mv.lor; 



end