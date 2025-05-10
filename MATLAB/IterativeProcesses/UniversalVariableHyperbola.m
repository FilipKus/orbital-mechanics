function [r,v] = UniversalVariableHyperbola(hmag,emag, TAinit, a, t, r0, v0)

% This function outputs the state vector (r,v) after time t given an
% initial state vector (r0, v0) for a given hyperbolic orbit.
%
% INPUTS:   hmag        - angular momentum
%           emag        - eccentricity
%           TAinit      - true anomaly (deg)
%           a           - semimajor-axis
%           t           - change in time
%           r0          - initial position vector of object (km)
%           v0          - initial velocity vector of object (km/s)
%
% OUTPUTS:  r           - final position vector of object (km)
%           v           - final velocity vector of object (km/s)
%
% Author: Filip Kus

mu = 398600;
vr0 = (mu/hmag)*emag*sin(TAinit);
rmag = norm(r0);
alpha = 1/a;

curX = sqrt(mu)* norm(alpha) * t;

condition = -1;

while condition == -1
    
    z = alpha * curX^2;
    Cz = (cosh(sqrt(-z)) - 1)/(-z);
    Sz = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z)^3);
    
    fX = (rmag*vr0/sqrt(mu)) * (curX^2) * Cz + (1 - alpha*rmag)*(curX^3) * Sz + rmag*curX - sqrt(mu)*t;    
    dfX = (rmag*vr0/sqrt(mu)) * curX *(1 - alpha*(curX^2)*Sz) + (1 - alpha*rmag)*(curX^2)*Cz + rmag;
    
    ratio = fX/dfX;
    
    curX = curX - ratio;
    if norm(ratio) < 10e-8
        condition = 1;
    end
end

f = 1 - (curX^2/rmag)*Cz;
g = t - (1/sqrt(mu)) * curX^3  * Sz;

r = f*r0 + g*v0;
rFinalMag = norm(r);

df = (sqrt(mu)/(rmag*rFinalMag)) * (alpha*curX^3 * Sz - curX);
dg = 1 - (curX^2 * Cz / rFinalMag);

v = df*r0 + dg*v0;
end        