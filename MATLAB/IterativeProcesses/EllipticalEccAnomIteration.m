function [Eone] = EllipticalEccAnomIteration(e, a, t, tp)

% This function outputs the eccentric anomaly used for an elliptical orbit.
%
% INPUTS:   e       - eccentricity
%           a       - semimajor-axis
%           t       - final time
%           tp      - time of perigee passage
%
% OUTPUTS:  Eone    - final iteration of the eccentric anomaly
%
% Author: Filip Kus

mu = 398600;
n = sqrt(mu/a^3);
M = n*(t-tp);
Ezero = M;

for i = 1:100
    Eone = Ezero - ((Ezero - e*sin(Ezero) - M)/(1 - e*cos(Ezero)));
    Ezero = Eone;
end
end


