clear;
clc;
% ----------------------------------------------------------------------%
% This script allows the user to convert between 3-2-1 and 3-1-3 euler
% angles to achieve the same rotations.
%
% Author: Filip Kus
% ----------------------------------------------------------------------%


% User selects conversion from 3-2-1 to 3-1-3 or the other way around
userChoice = input('Convert to 3-2-1 or 3-1-3: ','s');

%-----------------------------------
%-------Convert to 3-1-3 -----------
%-----------------------------------
if strcmp(userChoice, '3-1-3')
    
    % Preallocate array of angles
    angleSetDeg = zeros(1,3);    
    % User enters values of angles    
    angleSetDeg(1) = input('Enter the first angle (degrees): ');
    angleSetDeg(2) = input('Enter the second angle (degrees): ');
    angleSetDeg(3) = input('Enter the third angle (degrees): ');    
    % Convert to radians
    angleSetRad = angleSetDeg .* pi / 180;
    
    % Single-axis DCM's with each angle
    angleThreeDCM = [1, 0, 0;
               0, cos(angleSetRad(3)), sin(angleSetRad(3));
               0, -sin(angleSetRad(3)), cos(angleSetRad(3))];
           
    angleTwoDCM = [cos(angleSetRad(2)), 0, -sin(angleSetRad(2));
               0, 1, 0;
               sin(angleSetRad(2)), 0, cos(angleSetRad(2))];
         
    angleOneDCM = [cos(angleSetRad(1)), sin(angleSetRad(1)), 0;
               -sin(angleSetRad(1)), cos(angleSetRad(1)), 0;
               0, 0, 1];
    % Matrix multiplication to give DCM
    threeAxisDCM = angleThreeDCM * angleTwoDCM * angleOneDCM;
    
    % Calculations for each angle
    angleOne = rad2deg(atan(threeAxisDCM(3,1)/-threeAxisDCM(3,2)));
    angleTwo = rad2deg(acos(threeAxisDCM(3,3)));
    angleThree = rad2deg(atan(threeAxisDCM(1,3)/threeAxisDCM(2,3)));
    
    % Display results to user
    fprintf('The first angle is %.1f degrees\n', angleOne);
    fprintf('The second angle is %.1f degrees\n', angleTwo);
    fprintf('The third angle is %.1f degrees\n', angleThree);
end

%-----------------------------------
%-------Convert to 3-2-1 -----------
%-----------------------------------
if strcmp(userChoice, '3-2-1')
    
    % Preallocate array of angles
    angleSetDeg = zeros(1,3);    
    % User enters values of angles    
    angleSetDeg(1) = input('Enter the first angle (degrees): ');
    angleSetDeg(2) = input('Enter the second angle (degrees): ');
    angleSetDeg(3) = input('Enter the third angle (degrees): ');    
    % Convert to radians
    angleSetRad = angleSetDeg .* pi / 180;
    
    % Single-axis DCM's with each angle
    angleThreeDCM = [cos(angleSetRad(3)), sin(angleSetRad(3)), 0;
               -sin(angleSetRad(3)), cos(angleSetRad(3)), 0;
               0, 0, 1];
           
    angleTwoDCM = [1, 0, 0;
               0, cos(angleSetRad(2)), sin(angleSetRad(2));
               0, -sin(angleSetRad(2)), cos(angleSetRad(2))];
         
    angleOneDCM = [cos(angleSetRad(1)), sin(angleSetRad(1)), 0;
               -sin(angleSetRad(1)), cos(angleSetRad(1)), 0;
               0, 0, 1];
    % Matrix multiplication to give DCM
    threeAxisDCM = angleThreeDCM * angleTwoDCM * angleOneDCM;
    
    % Calculations for each angle
    angleOne = rad2deg(atan(threeAxisDCM(1,2)/threeAxisDCM(1,1)));
    angleTwo = rad2deg(asin(threeAxisDCM(1,3))*-1);
    angleThree = rad2deg(atan(threeAxisDCM(2,3)/threeAxisDCM(3,3)));
    
    % Display results to user
    fprintf('The first angle is %.1f degrees\n', angleOne);
    fprintf('The second angle is %.1f degrees\n', angleTwo);
    fprintf('The third angle is %.1f degrees\n', angleThree);
end

           
  