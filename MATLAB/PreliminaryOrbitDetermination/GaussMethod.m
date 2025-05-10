function [r2, v2] = GaussMethod(R1, R2, R3, rho1, rho2, rho3, t)

% This function outputs the state vectors r2 and v2 of an orbiting body 
% using Gauss's method with iterative improvement given three angle-only 
% observations at three closely-spaced times
%
% INPUTS:   R1      - observer position vector at time t1
%           R2      - observer position vector at time t2
%           R3      - observer position vector at time t3
%           rho1    - direction cosine vectors of satellite at t1
%           rho2    - direction cosine vectors of satellite at t2
%           rho3    - direction cosine vectors of satellite at t3
%           t       - 1x3 matrix containing times t1, t2 and t3
% 
% OUTPUTS:  r2      - position vector of satellite wrt geocentric
%                       equatorial at time t2
%           v2      - velocity vector of satellite wrt geocentric
%                       equatorial at time t2
%
% Author: Filip Kus


%-----------------------------------------------------------
%---------------- Gauss's method ---------------------------
%-----------------------------------------------------------
mu = 398600;

tInterval1 = t(1) - t(2);
tInterval3 = t(3) - t(2);
tInterval = tInterval3 - tInterval1;

p1 = cross(rho2, rho3);
p2 = cross(rho1, rho3);
p3 = cross(rho1, rho2);

D0 = dot(rho1, p1);

D = [dot(R1, p1) dot(R1, p2) dot(R1, p3);
    dot(R2, p1) dot(R2, p2) dot(R2, p3);
    dot(R3, p1) dot(R3, p2) dot(R3, p3)];

A = 1/D0*(-D(1,2)*tInterval3/tInterval + D(2,2) + D(3,2)*tInterval1/tInterval);

B = (1/(6*D0))*(D(1,2)*(tInterval3^2 - tInterval^2)*(tInterval3/tInterval) + D(3,2)*(tInterval^2 - tInterval1^2)*(tInterval1/tInterval));

E = dot(R2, rho2);
R2MagSq = dot(R2, R2);

alpha = -(A^2 + 2*A*E + R2MagSq);
b = -2*mu*B*(A + E);
c = -1*(mu^2)*B^2;

condition = -1;
x0 = 6800;
while condition == -1
    x1 = x0 - ((x0^8 + alpha*x0^6 + b*x0^3 + c)/(8*x0^7 + 6*alpha*x0^5 + 3*b*x0^2));
    
    if norm(x1 - x0) < 0.00001
        condition = 1;
    end   
    x0 = x1;
end

r2Mag = x1;
slantRange1 = 1/D0*((6*(D(3,1)*tInterval1/tInterval3 + D(2,1)*tInterval/tInterval3)*r2Mag^3 + mu*D(3,1)*(tInterval^2 - tInterval1^2)*tInterval1/tInterval3)/(6*r2Mag^3 + mu*(tInterval^2 - tInterval3^2)) - D(1,1));
slantRange2 = A + mu*B/r2Mag^3;
slantRange3 = 1/D0*((6*(D(1,3)*tInterval3/tInterval1 - D(2,3)*tInterval/tInterval1)*r2Mag^3 + mu*D(1,3)*(tInterval^2 - tInterval3^2)*tInterval3/tInterval1)/(6*r2Mag^3 + mu*(tInterval^2 - tInterval3^2)) - D(3,3));

r1 = R1 + slantRange1*rho1;
r2 = R2 + slantRange2*rho2;
r3 = R3 + slantRange3*rho3;

f1 = 1 - 0.5*mu*tInterval1^2/r2Mag^3;
f3 = 1 - 0.5*mu*tInterval3^2/r2Mag^3;
g1 = tInterval1 - (1/6)*mu*tInterval1^3/r2Mag^3;
g3 = tInterval3 - (1/6)*mu*tInterval3^3/r2Mag^3;

v2 = (-1*f3*r1 + f1*r3)/(f1*g3 - f3*g1);

% An approximate value for the state vectors r2 and v2 have been found. We
% now improve these values iteratively using Universal variables.

%-----------------------------------------------------------
%---------------- Iterative Improvement --------------------
%-----------------------------------------------------------

r2Mag = norm(r2);
r2MagF = 0;
iterations = 0;
while norm(r2MagF - r2Mag)> 0.1
    r2Mag = norm(r2);
    v2Mag = norm(v2);
    alpha = 2/r2Mag - v2Mag^2/mu;
    vRadial = dot(v2, r2)/r2Mag;

    curX1 = UniversalVarIteration(r2Mag, vRadial, tInterval1, alpha);
    curX3 = UniversalVarIteration(r2Mag, vRadial, tInterval3, alpha);

    z1 = alpha*curX1^2;
    f1I = 1- curX1^2/r2Mag*(1 - cos(sqrt(z1)))/z1;
    g1I = tInterval1 - curX1^3*(sqrt(z1) - sin(sqrt(z1)))/(sqrt(z1)^3)/sqrt(mu);
    z3 = alpha*curX3^2;
    f3I = 1- curX3^2/r2Mag*(1 - cos(sqrt(z3)))/z3;
    g3I = tInterval3 - curX3^3*(sqrt(z3) - sin(sqrt(z3)))/(sqrt(z3)^3)/sqrt(mu);

    f1 = (f1 + f1I)/2;
    f3 = (f3 + f3I)/2;
    g1 = (g1 + g1I)/2;
    g3 = (g3 + g3I)/2;

    c1 = g3/(f1*g3 - f3*g1);
    c3 = -1*g1/(f1*g3 - f3*g1);

    slantRange1 = (1/D0)*(-D(1,1) + D(2,1)/c1 - D(3,1)*c3/c1);
    slantRange2 = (1/D0)*(-c1*D(1,2) + D(2,2) - c3*D(3,2));
    slantRange3 = (1/D0)*(-D(3,3) + D(2,3)/c3 - D(1,3)*c1/c3);

    r1 = R1 + slantRange1*rho1;
    r2 = R2 + slantRange2*rho2;
    r3 = R3 + slantRange3*rho3;

    v2 = (-1*f3*r1 + f1*r3)/(f1*g3I - f3*g1);
    r2MagF = norm(r2);
    iterations = iterations + 1;
end
end



