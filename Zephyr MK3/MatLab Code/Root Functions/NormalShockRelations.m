function [MachOut, AxialPerUnitLength] = NormalShockRelations(MachIn,alpha)
%Mach MUST BE a column array Values greater than one
%Alpha MUST BE a row array IN DEGREES

% if you provide an incoming angle of attack array and mach array, the function will
%  provide an outgoing mach, outgoing pressure, and resulting axial force
% specifically from the normal shock flow and as an array or matrix.

% This assumes that the bow shock is normal and that PressureIn is always
% atmospheric.

thickness = 0.4*10^-2; % meters fixed for fins
gammaratio = 1.4; %Potentially a fixed value for air
PressureIn = 101325; % Pascals



%Rankine-Hugoniot Relations
    MachOut = sqrt((1+((gammaratio-1)/2).*MachIn.^2)./((gammaratio.*MachIn.^2)-((gammaratio-1)/2)));
    PressureOut = PressureIn * (1 + (2*gammaratio/(gammaratio + 1).*(MachIn.^2 - 1)));

%Force is pressure times area (thickness for 2d)
    AxialPerUnitLength = (PressureOut)*(thickness*secd(alpha)); %Matrix multiplication requires mach as column and alpha as row

%If want a coefficient instead, must divide by Dynamic Pressure
    
% AxialCoefficientPerUnitLength = AxialPerUnitLength/DynamicPressure;


end