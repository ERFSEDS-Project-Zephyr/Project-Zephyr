clear;
clc;
close all;
warning('off');

masterFile = readcell("ZephyrMK2_SimMasterFile.xlsx");
setup = cell2mat(masterFile(2:9,2));
bodyDims = cell2mat(masterFile(2:7,4));
bodyDims = cell2mat(masterFile(2:7+bodyDims(6),4));
finDims = cell2mat(masterFile(2:14,6:6+bodyDims(6)-1));
massDist = cell2mat(masterFile(2:25,10:12));

machMin = setup(1);
machMax = setup(2);
machDelta = setup(3);
attackAngMin = setup(4);
attackAngMax = setup(5);
attackAngDelta = setup(6);

bodyOuterDiameter = bodyDims(1)/100;
bodyInnerDiameter = bodyDims(2)/100;
noseLength = bodyDims(3)/100;
bodyDensity = bodyDims(4)*1000;
finResolution = bodyDims(5);
numStages = bodyDims(6);
bodyLength = sum(bodyDims(7:6+numStages))/100;

massDist = [massDist; [1000*bodyDensity*bodyLength*pi*0.25*(bodyOuterDiameter^2 - bodyInnerDiameter^2),100*bodyLength,100*(noseLength + bodyLength/2)]];
massDist(massDist(:,3)>100*(noseLength+bodyLength),:) = [];
massCG = massDist(:,1)/1000;
distCG = massDist(:,3)/100;
rocketCG = sum(massCG.*distCG)/sum(massCG);
distCG(end) = [];
massCG(end) = [];

k = 1;
for attackAng = (attackAngMin:attackAngDelta:attackAngMax)*pi/180
    j = 1;
    for Mach = machMin:machDelta:machMax
        fprintf('Mach: %.2f\t Attack Angle: %.2f deg.\n',Mach,attackAng*180/pi);
        sectionDist = [];
        sectionPointMass = [];
        if Mach <= 1
            [~,centerPres(:,j,k),liftForce(:,j,k),~] = subsonicRocket(bodyOuterDiameter,noseLength,numStages,Mach,attackAng,finDims,rocketCG);
        else
            [~,centerPres(:,j,k),liftForce(:,j,k),~] = supersonicRocket(bodyOuterDiameter,noseLength,numStages,finResolution,Mach,attackAng,finDims,rocketCG);
        end
        [sectionDist,sectionLengthIndex] = sort([distCG; centerPres(:,j,k)],'ascend');
        aeroIndex = find(sectionLengthIndex>length(distCG));
        for i = 1:length(sectionDist)
            if any(i == aeroIndex)
                sectionPointMass(i) = 0;
            else
                sectionPointMass(i) = massCG(sectionLengthIndex(i));
            end
        end
        sectionDist = [0;sectionDist;noseLength+bodyLength];
        sectionPointMass = [0;sectionPointMass';0];
        M = sysMassMat(sectionDist,sectionPointMass,bodyInnerDiameter,bodyOuterDiameter,bodyDensity);
        K = sysStiffMat(sectionDist,bodyInnerDiameter,bodyOuterDiameter,54*10^9,aeroIndex,liftForce(:,j,k),attackAng);
        D = M\K;
        modes(:,j,k) = sort(eig(D));
        if modes(2,j,k)<=0
            break;
        end
        safeMach(k,:) = [attackAng, Mach];
    
        j = j + 1;
    end
    clc;
    k = k + 1;
end

plot(safeMach);