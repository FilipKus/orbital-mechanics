function [Q_GH] = TopocentricGeoHor (latitude, theta)

% This function outputs the matrix of the transformation from geocentric
% equatorial (or topocentric equatorial) to the topocentric horizon
% coordiante system.
%
% INPUTS:   latitude        - latitude of the observer (deg)
%           theta           - local sidereal time of observer (deg)
%
% OUTPUTS:  Q_GH            - transformation matrix from geocentric
%                               equatorial to topocentric horizon.
%
% Author: Filip Kus

Q_HG = TopocentricHorGeo (latitude, theta);
Q_GH = transpose(Q_HG);

end