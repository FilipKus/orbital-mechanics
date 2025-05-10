function [r, v] = sv_from_coe(h, e, RAdeg, ideg, wdeg, TAdeg)

% This function outputs the state vector (r,v) from the classical orbital
% elements (coe).
%
% INPUTS:   h       - angular momentum
%           e       - eccentricity
%           RAdeg   - right ascension of the ascending node (deg)
%           ideg    - inclination of the orbit (deg)
%           wdeg    - argument of perigee (deg)
%           TAdeg   - true anomaly (deg)
%
% OUTPUTS:  r       - position vector in the geocentric equatorial frame (km)
%           v       - velocity vector in the geocentric equatorial frame (km/s)

mu = 398600;
RA      = deg2rad(RAdeg);
i       = deg2rad(ideg);
w       = deg2rad(wdeg);
TA      = deg2rad(TAdeg);




rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

R3_W = [cos(RA) sin(RA) 0;
        -sin(RA) cos(RA) 0;
            0       0    1];
        
R1_i = [1       0             0;
        0      cos(i)   sin(i);
        0     -sin(i)    cos(i)];

R3_w = [cos(w) sin(w) 0;
        -sin(w) cos(w) 0;
            0       0    1];
        
Q_pX = R3_W'*R1_i'*R3_w';

r = Q_pX * rp;
v = Q_pX * vp;

r = r';
v = v';
        
end

        
        
        
        
        
