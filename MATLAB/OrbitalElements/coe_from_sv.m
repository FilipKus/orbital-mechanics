function [hmag, emag, ideg, RAdeg, wdeg, TAdeg] = coe_from_sv(r, v)

% This function outputs the orbital elements given the state vectors radius
% and velocity.
%
% INPUTS:   r       - positional vector of object in orbit at time t
%           v       - velocity of object at time t
%
% OUTPUTS:  hmag    - angular momentum
%           emag    - eccentricty
%           ideg    - inclination
%           RAdeg   - right ascension of ascending node
%           wdeg    - argument of perigee
%           TAdeg   - true anomaly
%
% Author: Filip Kus

mu = 398600;
rmag = norm(r);

% Angular Momentum
h = cross(r, v);
hmag = norm(h);

% Eccentricity
e = cross(v,h)/mu - r/rmag;
emag = norm(e);

% Inclination
irad = acos(dot(h,[0,0,1])/hmag);
ideg = rad2deg(irad);

% Right ascension of the ascending node
N = cross([0,0,1], h);
Nmag = norm(N);
if Nmag == 0
    RAdeg = 0;
else
   RAtemp = acos(dot([1,0,0],N)/Nmag);
    if dot(N,[0,1,0]) >= 0
        RArad = RAtemp;
    else
        RArad = 2*pi - RAtemp;
    end

    RAdeg = rad2deg(RArad);
end


% Argument of Perigee
if Nmag == 0
    wdeg = 0;
    wrad = 0;
else
    wtemp = acos(dot(N,e)/(Nmag*emag));
    if dot(e, [0,0,1]) >= 0
        wrad = wtemp;
    else
        wrad = 2*pi - wtemp;
    end

    wdeg = rad2deg(wrad);
end


% True Anomaly
if Nmag == 0
    urad = 0;
else
    utemp = acos(dot(N,r)/(Nmag*rmag));
    if dot(r, [0,0,1]) >= 0
        urad = utemp;
    else
        urad = 2*pi - utemp;
    end
end

TArad = urad - wrad;
if TArad < 0 
    TArad = TArad + 2*pi;
end

TAdeg = rad2deg(TArad);
end
