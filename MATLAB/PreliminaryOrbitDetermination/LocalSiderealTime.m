function [theta] = LocalSiderealTime (d, m, y, hrs, mins, secs, eastLongitude)

% This function outputs the local sidereal time for a given date and time
% as well as location.
%
% INPUTS:   d               - current day
%           m               - current month
%           y               - current year
%           hrs             - current hour
%           mins           	- current minute
%           secs            - current second
%           eastLongitude   - angle between Greenwich and current location
%
% OUTPUTS:  theta           - local sidereal time (deg)
%
% Author: Filip Kus


J2000 = 2451545;
% Calculate Julian day number at 0 hr UT
J0 = 367*y - round(7*(y + round((m + 9)/12))/4) + round(275*m/9) + d + 1721013.5;

% Calculate time T0 in Julian centuries between the Julian day J0 and J2000
T0 = (J0 - J2000)/36525;

% Calculate the Greenwich sidereal time thetaG0 at 0 hr UT (deg)
thetaG0 = 100.4606184 + 36000.77004*T0 + 0.000387933*T0^2 - 2.583*10^-8 * T0^3;
thetaG0 = mod(thetaG0, 360);

% Calculate universal time (UT)
UT = hrs + mins/60 + secs/3600;

% Calculate the Greenwich sidereal time
thetaG = thetaG0 + 360.98564724*UT/24;
thetaG = mod(thetaG, 360);

% Calculate the local sidereal time
theta = thetaG + eastLongitude;
theta = mod(theta, 360);

end
