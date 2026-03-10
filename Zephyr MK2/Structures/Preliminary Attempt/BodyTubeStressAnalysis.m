clear; clc; close all; addpath('..'); addpath('../Aerodynamics');

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

maxAccel = 240; maxThrust = 3102;
numSectionsAxial = 1000; numSectionsRadial = 1000;
j = 1;
for Mach = machMin:machDelta:machMax
    k = 1;
    for attackAng = attackAngMin:attackAngDelta:attackAngMax
        if Mach==0
            normStress(:,:,1,:) = inf*ones(numSectionsAxial,numSectionsRadial,1,round((attackAngMax-attackAngMin)/attackAngDelta));
            shearStress(:,:,1,:) = inf*ones(numSectionsAxial,numSectionsRadial,1,round((attackAngMax-attackAngMin)/attackAngDelta));
            normFail(:,:,1,:) = zeros(numSectionsAxial,numSectionsRadial,1,round((attackAngMax-attackAngMin)/attackAngDelta));
            shearFail(:,:,1,:) = zeros(numSectionsAxial,numSectionsRadial,1,round((attackAngMax-attackAngMin)/attackAngDelta));
            break;
        elseif Mach<=1
            [~,centerPresX,liftForce,dragForce,~] = subsonicRocket(bodyOuterDiameter,noseLength,numStages,Mach,attackAng,finDims,rocketCG);
        else
            [~,centerPresX,liftForce,dragForce,~] = supersonicRocket(bodyOuterDiameter,noseLength,numStages,finResolution,Mach,attackAng,finDims,rocketCG);
        end
        dragForce(isnan(dragForce)) = 0;
        [normStress(:,:,j,k),shearStress(:,:,j,k),normFail(:,:,j,k),shearFail(:,:,j,k)] = internalStress(numSectionsAxial,numSectionsRadial,attackAng,bodyLength+centerPresX(1),bodyOuterDiameter,bodyInnerDiameter,maxAccel,maxThrust,distCG,massCG,centerPresX,liftForce,dragForce);

        if any(normFail,'all') || any(shearFail,'all')
            fprintf('Failure detected at M = %.2f, AOA = %.2f\n\n',Mach,attackAng*180/pi);
            break;
        end
        fprintf('M = %.2f\tAOA = %.2f\n',Mach,attackAng*180/pi);
        k = k + 1;
    end
    j = j + 1;
end