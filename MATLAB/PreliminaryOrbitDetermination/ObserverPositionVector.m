function [R] = ObserverPositionVector (latitude, theta, H) 

% This function outputs the position vector (as seen from the geocentric
% reference frame)of an observer given the latitude, local sidereal time 
% and elevation. 
%
% INPUTS:   latitude        - latitude of the observer (deg)
%           theta           - local sidereal time of observer (deg)
%           H               - elevation of observer (km)
%
% OUTPUTS:  R               - position vector of observer (km)
%
% Author: Filip Kus

Re = 6378;
f = 0.003353;

latitude = deg2rad(latitude);
theta = deg2rad(theta);

R = zeros(1,3);
R(1) = (Re/(sqrt(1 - (2*f - f^2)*(sin(latitude))^2)) + H)*cos(latitude)*cos(theta);
R(2) = (Re/(sqrt(1 - (2*f - f^2)*(sin(latitude))^2)) + H)*cos(latitude)*sin(theta);
R(3) = ((Re * (1-f)^2)/(sqrt(1 - (2*f - f^2)*(sin(latitude))^2)) + H)*sin(latitude);

end




