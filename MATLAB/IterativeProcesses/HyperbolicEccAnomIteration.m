function [Fone] = HyperbolicEccAnomIteration (e, h, t)

% This function outputs the hyperbolic eccentric anomaly.
%
% INPUTS:   e       - eccentricity
%           h       - angular momentum
%           t       - change in time
%
% OUTPUTS:  Fone    - final iteration of the eccentric anomaly
%
% Author: Filip Kus

mu = 398600;
Mh = t*(mu^2/h^3)*(e^2 - 1)^(3/2);
Fzero = Mh;

for i = 1:100
    Fone = Fzero - (e*sinh(Fzero) - Fzero - Mh)/(e*cosh(Fzero) - 1);
    Fzero = Fone;
end
end

