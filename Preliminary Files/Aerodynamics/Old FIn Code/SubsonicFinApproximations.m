%angle of attack on columns and mach number on row
%{
Using c_L = 2 * pi * alpha in degrees / sqrt(1-Mach^2)
Using c_D = (0.3 * t / chordRoot) + pi * (chordRoot + chordTip) / Span / (1-Mach^2) * angleAttack^2

s = height of one fin
A = area of the fin
t = thickness of the fin

AngleAttack = angle of attack in degree
minAngleAttack = starting Angle Attack from simMaster
maxAngleAttack = ending AngleAttack from simMaster
deltaAngleAttack = increment for the loop
Mach = speed in terms of mach
minMach = starting Mach value from simMaster
maxMach = ending Mach value must be less than 1
deltaMach = Mach increments
%}
%clear workspace, window, and allow access to outside functions
clear
clc
addpath('..')

%Create variables for the simMaster variables
simMaster = readcell("ZephyrMK2_SimMasterFile.xlsx");
setup = cell2mat(simMaster(2:9,2));
bodyDims = cell2mat(simMaster(2:7,4));
bodyDims = cell2mat(simMaster(2:7+bodyDims(6),4));
finDims = cell2mat(simMaster(2:15,6:6+bodyDims(6)-1));
massDist = cell2mat(simMaster(2:24,10:12));
%% Setup values from the simMaster
%Values unique to each fin set
chordRoot = finDims(2,:)/100; %convert from cm to m
chordTip = finDims(3,:)/100; %convert from cm to m
s = finDims(4,:)/100; %convert from cm to m
t = finDims(5,:)/100; %convert from cm to m

%Values same for all fin sets
minMach = setup(1);
deltaMach = setup(3); %changed from simMaster to be small
maxMach = setup(2)-deltaMach; %changed from simMaster to be below mach 1
minAngleAttack = setup(4); %leave in degrees
maxAngleAttack = setup(5);
deltaAngleAttack = setup(6);

%% Calculations for coefficients of lift and drag
%creates a 3 dimensional value where rows is mach, columns is angle of
%attack, and depth is finset
for finset = 1:bodyDims(6)
    j = 1;
    for AngleAttack = (minAngleAttack:deltaAngleAttack:maxAngleAttack)*pi/180 %converts to radians
        i = 1;
        for Mach = minMach:deltaMach:maxMach
            c_L(i,j,finset) = (2*pi*AngleAttack/sqrt(1-Mach^2));
            c_D(i,j,finset) = 0.3* t(finset) / chordRoot(finset) + pi * (chordRoot(finset) + chordTip(finset)) * (AngleAttack^2) / s(finset) / (1-Mach^2);
            i = i + 1;
        end
        j = j + 1;
    end
end


%Restate the coefficients per finset as a 2d matrix for cleaner viewing

c_Dfinset1 = c_D(:,:,1);
c_Dfinset2 = c_D(:,:,2);
%c_Dfinset3 = c_D(:,:,3); simMaster needs to state there are 3 stages

c_Lfinset1 = c_L(:,:,1);
c_Lfinset2 = c_L(:,:,2);
%c_Dfinset3 = c_L(:,:,3); simMaster needs to state there are 3 stages
%% Calculations for Forces and axial/normal coefficients
%{
Using LiftForce = q_inf * Span * (chordRoot + chordTip) * c_L (In Newtons)
Using DragForce = q_inf * Span * (chordRoot + chordTip) * c_D (In Newtons)

q_inf = shr / 2 * Pressure_inf in pascals * Mach^2
shr = specific heat constant pressure/ specific heat constant volume
      = 1.4 (a constant for air)
Pressure_inf = 101325 Pascals (the atmosphere)
%}

%creates a 3 dimensional value where rows is mach, columns is angle of
%attack, and depth is finset
shr = 1.4;
P_inf = 101325;
for finset = 1:bodyDims(6)
    j = 1;
    for AngleAttack = (minAngleAttack:deltaAngleAttack:maxAngleAttack)*pi/180
        i = 1;
        for Mach = minMach:deltaMach:maxMach
            q_inf = shr / 2 * P_inf * Mach^2;
            c_Lift(i,j,finset) = s(finset) * (chordRoot(finset) + chordTip(finset)) * c_L(i,j,finset) / chordRoot(finset);
            c_Drag(i,j,finset) = s(finset) * (chordRoot(finset) + chordTip(finset)) * c_D(i,j,finset) / chordRoot(finset);

            Ca_sub(i,j,finset) = c_Drag(i,j,finset)*cos(AngleAttack) + c_Lift(i,j,finset)*sin(AngleAttack);
            Cn_sub(i,j,finset) = -c_Drag(i,j,finset)*sin(AngleAttack) + c_Lift(i,j,finset)*cos(AngleAttack);
            i = i + 1;
        end
        j = j + 1;
    end
end


%Restate forces as 2d matrices for cleaner viewing (forces in Newtons)

Cafinset1 = Ca_sub(:,:,1);
Cafinset2 = Ca_sub(:,:,2);
% Cafinset3 = Ca(:,:,3); simMaster needs to say 3 stages

Cnfinset1 = Cn_sub(:,:,1);
Cnfinset2 = Cn_sub(:,:,2);
% Cnfinset3 = Cn(:,:,3); simMaster needs to say 3 stages

%% Saving variables to .mat

save("Subsonic_Axial",'Ca_sub')
save("Subsonic_Normal",'Cn_sub')
