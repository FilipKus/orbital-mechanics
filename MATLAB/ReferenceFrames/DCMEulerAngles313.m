clear;
clc;
% ----------------------------------------------------------------------%
% This script allows the user to convert between a set of 3-1-3 euler
% angles and a direction cosine matrix.
%
% Author: Filip Kus
% ----------------------------------------------------------------------%


% User inputs a DCM or three angles
userChoice = input('Find DCM or angles?: ','s');


%-----------------------------------
%-------Finding Angles--------------
%-----------------------------------
if strcmp(userChoice, 'angles')
    
    % Preallocate DCM
    DCM = zeros(3,3);
    % User enters values into DCM
    oneOne = input('Enter the (1,1) value: ');
    oneTwo = input('Enter the (1,2) value: ');
    oneThree = input('Enter the (1,3) value: ');
    twoThree = input('Enter the (2,3) value: ');
    threeThree = input('Enter the (3,3) value: ');
    
    % Calculations for each angle
    angleOne = rad2deg(atan(oneTwo/oneOne));
    angleTwo = rad2deg(asin(oneThree)*-1);
    angleThree = rad2deg(atan(twoThree/threeThree));
    
    % Display results to user
    fprintf('The first angle is %.1f degrees\n', angleOne);
    fprintf('The second angle is %.1f degrees\n', angleTwo);
    fprintf('The third angle is %.1f degrees\n', angleThree);
end


%-----------------------------------
%-------Finding DCM-----------------
%-----------------------------------
if strcmp(userChoice, 'DCM')
    
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
    twoAxisDCM = angleTwoDCM * angleOneDCM;
    singleAxisDCM = angleOneDCM;
    
    % Display result to user
    disp('The DCM is:')
    disp(threeAxisDCM);
    
    nCoordinate = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    
    firstBCoordinate = singleAxisDCM * nCoordinate;
    secondBCoordinate = twoAxisDCM * nCoordinate;
    thirdBCoordinate = threeAxisDCM * nCoordinate;
    
    % First rotation
    figure(1)
    hold on
    quiver3(0, 0, 0, 1, 0, 0, 'k')
    hold on
    quiver3(0, 0, 0, 0, 1, 0,'k')
    hold on
    quiver3(0, 0, 0, 0, 0, 1,'k')
    hold on
    
    quiver3(0, 0, 0, firstBCoordinate(1,1), firstBCoordinate(1,2), firstBCoordinate(1,3), 'r')
    hold on
    quiver3(0, 0, 0, firstBCoordinate(2,1), firstBCoordinate(2,2), firstBCoordinate(2,3),'r')
    hold on
    quiver3(0, 0, 0, firstBCoordinate(3,1), firstBCoordinate(3,2), firstBCoordinate(3,3),'r')
    hold off
    view(140,25)
    
    
    % Second rotation
    figure(2)
    hold on
    quiver3(0, 0, 0, 1, 0, 0, 'k')
    hold on
    quiver3(0, 0, 0, 0, 1, 0,'k')
    hold on
    quiver3(0, 0, 0, 0, 0, 1,'k')
    hold on
    
    quiver3(0, 0, 0, secondBCoordinate(1,1), secondBCoordinate(1,2), secondBCoordinate(1,3), 'r')
    hold on
    quiver3(0, 0, 0, secondBCoordinate(2,1), secondBCoordinate(2,2), secondBCoordinate(2,3),'r')
    hold on
    quiver3(0, 0, 0, secondBCoordinate(3,1), secondBCoordinate(3,2), secondBCoordinate(3,3),'r')
    hold off
    view(140,25)
    
        
    % Third rotation
    figure(3)
    quiver3(0, 0, 0, 1, 0, 0, 'k')
    hold on
    quiver3(0, 0, 0, 0, 1, 0,'k')
    hold on
    quiver3(0, 0, 0, 0, 0, 1,'k')
    hold on

    quiver3(0, 0, 0, thirdBCoordinate(1,1), thirdBCoordinate(1,2), thirdBCoordinate(1,3), 'r')
    hold on
    quiver3(0, 0, 0, thirdBCoordinate(2,1), thirdBCoordinate(2,2), thirdBCoordinate(2,3),'r')
    hold on
    quiver3(0, 0, 0, thirdBCoordinate(3,1), thirdBCoordinate(3,2), thirdBCoordinate(3,3),'r')
    hold off
    view(140,25)
end




          
