function [azimuthDeg, angularElevationDeg] = AzimuthAngularElevation (sigma)

% This function outputs the azimuth and angular elevation relative to an
% observer.
%
% INPUTS:   sigma                   - position vector of object relative to
%                                       observer (km)
%
% OUTPUTS:  azimuthDeg              - topocentric right ascension (deg)
%           angularElevationDeg     - topocentric declination (deg)
%
% Author: Filip Kus


% Obtain unit vectors from the position of the object
sigmaUnitVec = sigma/norm(sigma);

% Obtain the topocentric declension
angularElevation = asin(sigmaUnitVec(3));

% Obtain the right ascension
azimuth = acos(sigmaUnitVec(2)/cos(angularElevation));

if sigmaUnitVec(1)/cos(angularElevation) < 0
    azimuth = 2*pi - azimuth;
end

% Convert to degrees
angularElevationDeg = rad2deg(angularElevation);
azimuthDeg = rad2deg(azimuth);

end
