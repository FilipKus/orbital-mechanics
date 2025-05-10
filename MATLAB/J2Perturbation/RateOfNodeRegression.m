function [nodeRegressionRad, nodeRegressionDeg] = RateOfNodeRegression (J2, R, e, a, ideg)

% This function calculates the rate of node regression.
% INPUTS: J2        - dimensionless parameter quantifying effects of oblatness of
%                       planetary body on orbit.
%         R         - radius of body
%         e         - eccentricity
%         a         - semimajor-axis
%         ideg      - inclination of orbit (deg)
%
% OUTPUTS: nodeRegressionRad - rate of node regression (rad/s)                             
%          nodeRegressionDeg - rate of node regression (deg/day)
%
% Author: Filip Kus

mu = 398600;
irad = deg2rad(ideg);

nodeRegressionRad = -1*cos(irad)*(3/2)*sqrt(mu)*J2*R^2 / ((1-e^2)^2 * a^(7/2));
nodeRegressionDeg = nodeRegressionRad * (180/pi) * 86400;
end
