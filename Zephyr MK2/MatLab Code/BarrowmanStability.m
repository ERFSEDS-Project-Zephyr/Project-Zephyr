clear; clc; close all; warning('off','all'); addpath('..'); addpath('Barrowman Functions')
%<Read and parse info in file>
simMaster = readcell("ZephyrMK2_SimMasterFile.xlsx");
setup = cell2mat(simMaster(2:11,2));
bodyDims = cell2mat(simMaster(2:7,4));
bodyDims = cell2mat(simMaster(2:7+bodyDims(6),4));
finDims = cell2mat(simMaster(3:16,6:6+bodyDims(6)-1));
massDist = cell2mat(simMaster(3:18,9:11));

machMin = setup(3);
machMax = setup(4);
machDelta = setup(5);
attackAngMin = setup(6)*pi/180;
attackAngMax = setup(7)*pi/180;
attackAngDelta = setup(8)*pi/180;
stabilityMin = setup(9);
stabilityMax = setup(10);

bodyOuterDiameter = bodyDims(1)/100;
bodyInnerDiameter = bodyDims(2)/100;
noseLength = bodyDims(3)/100;
bodyDensity = bodyDims(4)*1000;
finResolution = bodyDims(5);
numStages = bodyDims(6);
bodyLength = sum(bodyDims(7:6+numStages))/100;

massDist = [massDist; [1000*bodyDensity*bodyLength*pi*(bodyOuterDiameter^2 - bodyInnerDiameter^2)/4,100*bodyLength,100*(noseLength + bodyLength/2)]];
massDist(massDist(:,3)>100*(noseLength+bodyLength),:) = [];
massCG = massDist(:,1)/1000;
distCG = massDist(:,3)/100;
rocketCG = sum(distCG.*massCG)/sum(massCG);
%<\Read and parse info in file>

%<Compute stability>
j=1;
for Mach = machMin:machDelta:machMax
    k=1;
    for attackAng = attackAngMin:attackAngDelta:attackAngMax
        if Mach<=1
            [~,~,liftForce(k,j,:),dragForce(k,j,:),stability(k,j)] = subsonicRocket(bodyOuterDiameter,noseLength,bodyLength,numStages,Mach,attackAng,finDims,rocketCG);
        else
            [~,~,liftForce(k,j,:),dragForce(k,j,:),stability(k,j)] = supersonicRocket(bodyOuterDiameter,noseLength,bodyLength,numStages,finResolution,Mach,attackAng,finDims,rocketCG);
        end
        k=k+1;
    end
    fprintf("\b\b\b\b\b\b\b\b%.2f%%",100*j/(1 + (machMax-machMin)/machDelta));
    j=j+1;
end
fprintf("\n");
stability(stability<stabilityMin) = stabilityMin;
stability(stability>stabilityMax) = stabilityMax;
%</Compute stability>

attackAng = (attackAngMin:attackAngDelta:attackAngMax)*180/pi;
Mach = machMin:machDelta:machMax;

%<2D plot stability near 0 deg. angle of attack>
title("Stability and CP vs Mach @ ~0 Angle of Attack");
plot(Mach,stability(2,:),'LineWidth',3);
xlabel("Mach Number");
ylabel("Stability (cal)");
xlim([machMin,machMax]);
ylim([stabilityMin,stabilityMax]);
%</2D plot stability near 0 deg. angle of attack>