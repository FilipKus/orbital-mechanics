function[r,v] = PerifocalToGeocentric(r0,v0,ideg,RAdeg,wdeg)

% This function converts the position and velocity vectors of an object in
% orbit from a perifocal frame of reference to a geocentric frame of
% reference.
%
% INPUTS:   r0      - position vector in a perifocal reference frame (km)
%           v0      - velocity vector in a perifocal reference frame (km/s)
%           ideg    - inclination (deg)
%           RAdeg   - right ascencion of ascending node (deg)
%           wdeg    - argument of perigee (deg)
%
% OUTPUTS:  r       - position vector in a geocentric reference frame (km)
%           v       - velocity vector in a geocentric reference frame (km/s)
%
% Author: Filip Kus

i = deg2rad(ideg);
RA = deg2rad(RAdeg);
w = deg2rad(wdeg);

Q = [(cos(w)*cos(RA)-sin(w)*cos(i)*sin(RA))   (cos(w)*sin(RA)+sin(w)*cos(i)*cos(RA))     sin(w)*sin(i);
    (-sin(w)*cos(RA)-cos(w)*cos(i)*sin(RA))     (-sin(w)*sin(RA)+cos(w)*cos(i)*cos(RA))     cos(w)*sin(i);
    sin(i)*sin(RA)                                            -sin(i)*cos(RA)               cos(i)];


Qprime = transpose(Q);

r = Qprime * transpose(r0);
v = Qprime * transpose(v0);
end
