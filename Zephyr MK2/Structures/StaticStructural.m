clear; clc; close all; addpath('..'); addpath('../Aerodynamics');

%<Read and parse info in file>
simMaster = readcell("ZephyrMK2_SimMasterFile.xlsx");
setup = cell2mat(simMaster(2:13,2));
bodyDims = cell2mat(simMaster(2:7,4));
bodyDims = cell2mat(simMaster(2:7+bodyDims(6),4));
finDims = cell2mat(simMaster(3:16,6:6+bodyDims(6)-1));
massDist = cell2mat(simMaster(3:18,9:11));

gamma = setup(1);
R = setup(2);
pinf = setup(3)*10^5;
Tinf = setup(4) + 273.15;
machMin = setup(5);
machMax = setup(6);
machDelta = setup(7);
attackAngMin = setup(8)*pi/180;
attackAngMax = setup(9)*pi/180;
attackAngDelta = setup(10)*pi/180;
stabilityMin = setup(11);
stabilityMax = setup(12);

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

I = (pi/64)*(bodyOuterDiameter^4 - bodyInnerDiameter^4);
J = (pi/32)*(bodyOuterDiameter^4 - bodyInnerDiameter^4);
area_cs = (pi/4)*(bodyOuterDiameter^2 - bodyInnerDiameter^2);

m = 1;
for M = machMin:machDelta:machMax    
    a = 1;
    for alpha = attackAngMin:attackAngDelta:attackAngMax
        [A(m,a,1),N(m,a,1),xcp(m,a,1),ycp(m,a,1)] = noseAero(gamma,R,pinf,M,alpha,bodyOuterDiameter,noseLength,0.5*(bodyOuterDiameter-bodyInnerDiameter));
        for n = 1:numStages
            [A(m,a,n+1),N(m,a,n+1),xcp(m,a,n+1),ycp(m,a,n+1)] = finAero(gamma,R,pinf,Tinf,M,alpha,bodyOuterDiameter,finDims(:,n),finResolution);
        end
        
        
        
        a = a + 1;
    end
    m = m + 1;
end

