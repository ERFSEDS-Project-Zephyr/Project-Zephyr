%Clear workspace
clear;
close all;
clc;
addpath('../Data Files');

%Read data from sim file
simMaster = readcell("ZephyrMK3_SimMasterFile.xlsx");
setup = cell2mat(simMaster(2:9,2));
bodyDims = cell2mat(simMaster(2:8,4));
bodyDims = cell2mat(simMaster(2:8+bodyDims(7),4));
finDims = cell2mat(simMaster(3:11,6:6+bodyDims(7)-1));


%Values unique to each fin set
chordRoot = finDims(2,:)/100; %convert from cm to m
chordTip = finDims(3,:)/100; %convert from cm to m
s = finDims(4,:)/100; %convert from cm to m
t = finDims(5,:)/100; %convert from cm to m



%Hardcoded values
minMach=0;
gamma=1.4;
minAngleAttack=0;
maxAngleAttack=5;
DeltaAngleAttack=0.1;
maxMach=1;
deltaMach=0.01
p_inf=101325;

MachVect = minMach:deltaMach:maxMach;
alphaVect = minAngleAttack:DeltaAngleAttack:maxAngleAttack;
NormForceTable = zeros(length(MachVect),length(alphaVect)); 
AxialForceTable = zeros(length(MachVect),length(alphaVect)); 

%Calculations

for finset = 1:bodyDims(7)

    j = 1;

    for AngleAttack = deg2rad((minAngleAttack:DeltaAngleAttack:maxAngleAttack))
        
        i = 1;
       
        for Mach = minMach:deltaMach:maxMach-deltaMach
            if Mach <= 1-deltaMach
                c_bar(1,1:finset)=(chordTip(finset)+chordRoot(finset))/2;
                c_Di(1,1:finset)=(2*pi*(c_bar/s(finset)))*(AngleAttack.^2)/(1+(c_bar/s(finset)).^2);
                C_mac(1,1:finset)=(4/(3.*c_bar(1,finset))).*((c_bar(1,finset).^2)-(chordTip(finset).*chordRoot(finset)));
                c_D = c_Di + (0.3 .* (t(finset)/ C_mac(1,finset)).^2); 
                c_L=(2*pi*AngleAttack)/(1+(c_bar(1,finset)/s(finset)));
                q_inf=(gamma/2)*p_inf*(Mach)^2;

                D=c_D(1,finset)*q_inf*s(1,finset)*c_bar(1,finset);
                L=c_L*q_inf*s(1,finset)*c_bar(1,finset);
                NormForce = L*cos(AngleAttack)+D*sin(AngleAttack);
                Axial = D*cos(AngleAttack)+L*sin(AngleAttack);
        
                
                Cn(i,j,1:finset) = -c_D.*sin(AngleAttack) + c_L.*cos(AngleAttack);
                Ca(i,j,1:finset) = (c_D.*cos(AngleAttack) + c_L.*sin(AngleAttack))/sqrt((1-(Mach^2)));
            end
            i=i+1;
        end
        j = j + 1; % Increment angle attack index
    end
end
