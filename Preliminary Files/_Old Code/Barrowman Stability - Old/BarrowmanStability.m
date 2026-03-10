clear;
clc;
close all;
warning('off','all');
%<Read and parse info in file>
simMaster = readcell("ZephyrMK2_SimMasterFile.xlsx");
setup = cell2mat(simMaster(2:9,2));
bodyDims = cell2mat(simMaster(2:7,4));
bodyDims = cell2mat(simMaster(2:7+bodyDims(6),4));
finDims = cell2mat(simMaster(2:15,6:6+bodyDims(6)-1));
massDist = cell2mat(simMaster(2:24,10:12));

machMin = setup(1);
machMax = setup(2);
machDelta = setup(3);
attackAngMin = setup(4)*pi/180;
attackAngMax = setup(5)*pi/180;
attackAngDelta = setup(6)*pi/180;
stabilityMin = setup(7);
stabilityMax = setup(8);

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
    for attackAng=attackAngMin:attackAngDelta:attackAngMax
        if Mach<=1
            [~,~,~,~,stability(k,j)] = subsonicRocket(bodyOuterDiameter,noseLength,numStages,Mach,attackAng,finDims,rocketCG);
        else
            [~,~,~,~,stability(k,j)] = supersonicRocket(bodyOuterDiameter,noseLength,numStages,finResolution,Mach,attackAng,finDims,rocketCG);
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