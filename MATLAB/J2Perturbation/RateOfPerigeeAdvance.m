function [perigeeAdvanceRad, perigeeAdvanceDeg] = RateOfPerigeeAdvance (J2, R, e, a, ideg)

% This function calculates the rate of perigee advance.
% INPUTS:   J2      - dimensionless parameter quantifying effects of oblatness of
%                       planetary body on orbit.
%           R       - radius of body
%           e       - eccentricity
%           a       - semimajor-axis
%           ideg    - inclination of orbit (deg)
%
% OUTPUTS:  perigeeAdvanceRad - rate of perigee advance (rad/s)
%           perigeeAdvanceDeg - rate of perigee advance (deg/day)
%
% Author: Filip Kus

mu = 398600;
irad = deg2rad(ideg);

perigeeAdvanceRad = -1*((5/2)*(sin(irad))^2 - 2)*(3/2)*sqrt(mu)*J2*R^2 / ((1-e^2)^2 * a^(7/2));
perigeeAdvanceDeg = perigeeAdvanceRad * (180/pi) * 86400;
end

