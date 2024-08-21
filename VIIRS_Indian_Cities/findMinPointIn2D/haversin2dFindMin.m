% MATLAB code for the minimum distance approximated by Haversine Formula, which can be found at the link below: 
% https://en.wikipedia.org/wiki/Haversine_formula
function [I, J, haversine_ch] = haversin2dFindMin(LatInDecDegrees, LonInDecDegrees, latInDecDegrees, lonInDecDegrees, sphereRadiusInDesiredUnits)
    phi1 = (pi/180) * latInDecDegrees;
    phi2 = (pi/180) * LatInDecDegrees;
    delta_lat=deg2rad(LatInDecDegrees - latInDecDegrees);   
    delta_lon=deg2rad(LonInDecDegrees - lonInDecDegrees);                   
    a = sin((delta_lat./2).^2) + cos(phi1) .* cos(phi2) .* sin((delta_lon./2).^2);
    a(a < 0) = 0; a(a > 1) = 1;    
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    distances = sphereRadiusInDesiredUnits * c;
    [I,J]=find(distances(:,:)==min(distances(:)));
    haversine_ch = distances(I,J);  	
end