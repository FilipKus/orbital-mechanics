function [a, e, ideg, RAdeg, wdeg, TAdeg] = GibbsMethod(r1, r2, r3)

% This function uses Gibb's method of preliminary orbit determination to
% determine the classical orbital elements for a given set of three
% geocentric position vectors of a space object at three successive times.
%
% INPUTS:   r1      - geocentric position vector at time t1
%           r2      - geocentric position vector at time t2
%           r3      - geocentric position vector at time t3
%
% OUTPUTS:  a       - semi-major axis 
%           e       - eccentricty
%           ideg    - inclination (deg)
%           RAdeg   - right ascension of the ascending node (deg)
%           wdeg    - argument of perigee (deg)
%           TAdeg   - true anomaly at time t2(deg)
%
% Author: Filip Kus

mu = 398600;
%---------------------------------------------------%
%-------------- Finding v2 -------------------------%
%---------------------------------------------------%
% Find the magnitude of the position vectors
r1mag = norm(r1);
r2mag = norm(r2);
r3mag = norm(r3);

% Find the coefficients
C12 = cross(r1, r2);
C23 = cross(r2, r3);
C31 = cross(r3, r1);

% Calculated N, D and S
N = r1mag*C23  + r2mag*C31 + r3mag*C12;
D = C12 + C23 + C31;
S = r1*(r2mag - r3mag) + r2*(r3mag - r1mag) + r3*(r1mag - r2mag);

Nmag = norm(N);
Dmag = norm(D);

% Find v2
v2 = sqrt(mu/(Nmag*Dmag))*((cross(D, r2)/r2mag) + S);

%---------------------------------------------------%
%-------------- Orbital Elements--------------------%
%---------------------------------------------------%

[hmag, emag, ideg, RAdeg, wdeg, TAdeg] = coe_from_sv(r2, v2);

a = (hmag^2/mu)/(1-emag^2);
e = emag;
end



