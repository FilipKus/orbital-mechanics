function [curX] = UniversalVarIteration (r2Mag, vRadial, dt, alpha)
% This function outputs the universal anomaly by using Newton's method to
% solve the universal Kepler equation.
%
% INPUTS:   r2Mag   - position of satellite
%           vRadial - radial velocity
%           dt      - change in time
%           alpha   - reciprocal of the semimajor-axis
% 
% OUTPUTS:  curX    - universal anomaly
%
% Author: Filip Kus

mu = 398600;
curX = sqrt(mu)* norm(alpha) * dt;
condition = -1;
while condition == -1
    
    z = alpha * curX^2;
    Cz = (1 - cos(sqrt(z)))/z;
    Sz = (sqrt(z) - sin(sqrt(z)))/(sqrt(z)^3);
    
    fX = (r2Mag*vRadial/sqrt(mu)) * (curX^2) * Cz + (1 - alpha*r2Mag)*(curX^3) * Sz + r2Mag*curX - sqrt(mu)*dt;    
    dfX = (r2Mag*vRadial/sqrt(mu)) * curX *(1 - alpha*(curX^2)*Sz) + (1 - alpha*r2Mag)*(curX^2)*Cz + r2Mag;
    
    ratio = fX/dfX;
    
    curX = curX - ratio;
    if norm(ratio) < 10e-8
        condition = 1;
    end
end
end
