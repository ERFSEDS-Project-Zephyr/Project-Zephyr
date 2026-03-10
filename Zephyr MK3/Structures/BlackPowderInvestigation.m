clear; close all; clc; warning('off'); addpath('../..');
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',24,'defaultAxesFontName','Times New Roman');

% Givens
outerDiameter = 3.12*25.4; % Body Outer Diameter (mm)
innerDiameter = 3.00*25.4; % Body Inner Diameter (mm)
delta = outerDiameter/innerDiameter;
yieldX = 125; % Body Tube Longitudinal Yield Strength (MPa)
yieldY = 45.4; % Body Tube Tangential Yield Strength (MPa)

% Simulation Parameters
numPins = 8; % Number of Shear Pins
pinDiameter = convert(0.1380,'in','mm'); % Shear Pin Diameter (mm)
yieldPins = convert(114,'lbf','N')/(0.25*pi*(pinDiameter^2)); % Shear Pin Yield Strength (MPa)

alpha = (yieldX-yieldY)*(yieldX^-2 - yieldY^-2);
yieldZ = (1 + sign(alpha)*sqrt(1 - alpha*(yieldX + yieldY)))/alpha;
H = 0.5*(yieldX^-2 + yieldY^-2 - yieldZ^-2);
p_BT = (delta^2 - 1)*(delta - 1)/sqrt(((delta-1)/yieldX)^2 + ((delta^2 - 1)/yieldY)^2 - 2*H*(delta^2 - 1)*(delta - 1));
p_SP = 2*yieldPins/sqrt((delta/(delta^2 - 1))^2 + (4/(numPins^2))*((innerDiameter/pinDiameter)^4));

if p_SP < p_BT
    fprintf('Body tube will NOT fail (F.S. = %.2f)\n',p_BT/p_SP);
else
    fprintf('Body tube will fail (F.S. = %.2f)\n',p_BT/p_SP);
end