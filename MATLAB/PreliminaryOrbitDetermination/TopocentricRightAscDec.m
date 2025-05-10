function [rightAscensionDeg, declinationDeg] = TopocentricRightAscDec (sigma)

% This function outputs the topocentric right ascension and declination of
% an observer given the latitude, local sidereal time and elevation. 
%
% INPUTS:   sigma               - position vector of object relative to
%                                   observer (km)
%
% OUTPUTS:  rightAscensionDeg   - topocentric right ascension (deg)
%           declinationDeg      - topocentric declination (deg)
%
% Author: Filip Kus


% Obtain unit vectors from the position of the object
sigmaUnitVec = sigma/norm(sigma);

% Obtain the topocentric declension
declination = asin(sigmaUnitVec(3));

% Obtain the right ascension
rightAscension = acos(sigmaUnitVec(1)/cos(declination));

if sigmaUnitVec(2)/cos(declination) < 0
    rightAscension = 2*pi - rightAscension;
end

% Convert to degrees
declinationDeg = rad2deg(declination);
rightAscensionDeg = rad2deg(rightAscension);

end
