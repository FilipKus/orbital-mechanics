function [a, e, ideg, RAdeg, wdeg, TAdeg] = LambertsProblem(r1, r2, t, trajectory)

% This function tackles Lambert's problem using universal variables given
% two position vectors, time between these two vectors and the trajectory
% of the orbit.
%
% INPUTS:   r1      - position vector at time t1
%           r2      - position vector at time t2
%           t       - time elapsed between the two position vectors.
%           trajectory - the inclination of the orbit (prograde or
%                           retrograde)
%
% OUTPUTS:  a       - semi-major axis 
%           e       - eccentricty
%           ideg    - inclination (deg)
%           RAdeg   - right ascension of the ascending node (deg)
%           wdeg    - argument of perigee (deg)
%           TAdeg   - true anomaly at time t1 (deg)
%
% Author: Filip Kus

mu = 398600;

% Obtain magnitude of position vectors
r1mag = norm(r1);
r2mag = norm(r2);

% Determine change in true anomaly based on the z component of the cross
% product of the two position vectors and whether the orbit is prograde or
% retrograde.

crossr12 = cross(r1, r2);
if trajectory == "prograde"
    if crossr12(3) >= 0
        deltaTA = acos(dot(r1,r2)/(r1mag*r2mag));        
    else
        deltaTA = 2*pi - acos(dot(r1,r2)/(r1mag*r2mag));        
    end
end
if trajectory == "retrograde"
    if crossr12(3) < 0
        deltaTA = acos(dot(r1,r2)/(r1mag*r2mag));        
    else
        deltaTA = 2*pi - acos(dot(r1,r2)/(r1mag*r2mag));        
    end
end

A = sin(deltaTA)*sqrt((r1mag*r2mag)/(1-cos(deltaTA)));

% By iteration using Newton's method, determine the value for z and thus
% determine whether the orbit is elliptical (z > 0) or hyperbolic (z < 0)
condition = -1;
z = 0;
while condition == -1
    if z > 0        
        Cz = (1 - cos(sqrt(z)))/z;
        Sz = (sqrt(z) - sin(sqrt(z)))/(sqrt(z)^3);
        yz = r1mag + r2mag + A * ((z*Sz - 1)/sqrt(Cz));        
        Fz = Sz * (yz/Cz)^(3/2) + A*sqrt(yz) - t*sqrt(mu);
        dFz = (sqrt(2)/40)*yz^(3/2) + (A/8)*(sqrt(yz) + A*sqrt(1/(2*yz)));
    
        ratio = Fz/dFz;
        
        z = z - ratio;
    elseif z < 0                
        Cz = (cosh(sqrt(-z)) - 1)/(-z);
        Sz = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z)^3);
        yz = r1mag + r2mag + A * ((z*Sz - 1)/sqrt(Cz));        
        Fz = Sz * (yz/Cz)^(3/2) + A*sqrt(yz) - t*sqrt(mu);
        dFz = (sqrt(2)/40)*yz^(3/2) + (A/8)*(sqrt(yz) + A*sqrt(1/(2*yz)));
    
        ratio = Fz/dFz;
        
        z = z - ratio;
    elseif z == 0       
        Cz = 1/2;
        Sz = 1/6;
        yz = r1mag + r2mag + A * ((z*Sz - 1)/sqrt(Cz));
        Fz = Sz * (yz/Cz)^(3/2) + A*sqrt(yz) - t*sqrt(mu);
        dFz = (sqrt(2)/40)*yz^(3/2) + (A/8)*(sqrt(yz) + A*sqrt(1/(2*yz)));
        
        ratio = Fz/dFz;
        
        z = z - ratio;
    end      
    
    if norm(ratio) < 10e-8
        condition = 1;
    end
end

% Calculate the Lagrange functions f, g and dg. 
f = 1 - yz/r1mag;
g = A * sqrt(yz/mu);
dg = 1 - yz/r2mag;

% Calculate the velocity vectors at time t1 and t2
v1 = (1/g)*(r2 - f*r1);
v2 = (1/g)*(dg*r2 - r1);

% Obtain the classical orbital elements
[hmag, emag, ideg, RAdeg, wdeg, TAdeg] = coe_from_sv(r1, v1);
a = (hmag^2/mu)/(1-emag^2);
e = emag;
end

