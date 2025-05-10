function [Q_HG] = TopocentricHorGeo (latitude, theta)

% This function outputs the matrix of the transformation from topocentric
% horizon to geocentric equatorial coordinate sytem (or topocentric
% equatorial).
%
% INPUTS:   latitude        - latitude of the observer (deg)
%           theta           - local sidereal time of observer (deg)
%
% OUTPUTS:  Q_HG            - transformation matrix from topocentric
%                               horizon to geocentric equatorial.
%
% Author: Filip Kus

latitude = deg2rad(latitude);
theta = deg2rad(theta);

Q_HG = [-sin(theta) -sin(latitude)*cos(theta) cos(latitude)*cos(theta);
        cos(theta)  -sin(latitude)*sin(theta)   cos(latitude)*sin(theta);
        0           cos(latitude)               sin(latitude)];

end