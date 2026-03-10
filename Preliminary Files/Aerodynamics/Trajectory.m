clear;
clc;
close all;
addpath('..');

%<Extract rocket geometry and thrust-mass data>
simMaster = readcell("ZephyrMK2_SimMasterFile.xlsx");
motorThrustMass = uigetfile("*.xlsx");
motorThrustMass = readcell(motorThrustMass);

setup = cell2mat(simMaster(2:13,2));
bodyDims = cell2mat(simMaster(2:7,4));
bodyDims = cell2mat(simMaster(2:7+bodyDims(6),4));
finDims = cell2mat(simMaster(3:16,6:5+bodyDims(6)));
massDist = cell2mat(simMaster(3:18,9:11));

machMin = setup(5);
machMax = setup(6);
machDelta = setup(7);
attackAngMin = setup(8)*pi/180;
attackAngMax = setup(9)*pi/180;
attackAngDelta = setup(10)*pi/180;

bodyOuterDiameter = bodyDims(1)/100;
bodyInnerDiameter = bodyDims(2)/100;
noseLength = bodyDims(3)/100;
bodyDensity = bodyDims(4)*1000;
finResolution = bodyDims(5);
numStages = bodyDims(6);
bodyLength = sum(bodyDims(6:5+numStages))/100;

massDist = [massDist; [1000*bodyDensity*bodyLength*pi*(bodyOuterDiameter^2 - bodyInnerDiameter^2)/4,100*bodyLength,100*(noseLength + bodyLength/2)]];
massDist(massDist(:,3)>100*(noseLength+bodyLength),:) = [];
massCG = massDist(:,1)/1000;
distCG = massDist(:,3)/100;
rocketCG = sum(distCG.*massCG)/sum(massCG);

timeStamp = cell2mat(motorThrustMass(2:end,1));
Thrust = cell2mat(motorThrustMass(2:end,2));
Mass = cell2mat(motorThrustMass(2:end,3));
MassDot = cell2mat(motorThrustMass(2:end,4));
%</Extract rocket geometry and thrust-mass data>

%<Establish initial conditions and universal constants>
T0 = 304.8;
T = T0;
rho0 = 1.225;
rho = rho0;
h = 0;
v = 0;
a = 0;
Mach = 0;
[~,~,~,~,stability] = subsonicRocket(bodyOuterDiameter,noseLength,bodyLength,2,Mach,0.01,finDims,rocketCG);
Cd = 0;
dt = 0.001;
t = dt;

g = 9.807;
B = 0.0065;
gamma = 1.4;
R = 8.314/0.0288;
%</Establish initial conditions and universal constants>

numFinSets = 2;

while t<=max(timeStamp)||v(end)>0
    %<Stage separation>
    if any(t == round(timeStamp(MassDot < -5),3))
        numStages = numStages - 1;
        bodyDims(end) = [];
        bodyLength = sum(bodyDims(7:6+numStages))/100;
        rocketLength = noseLength + bodyLength;
        massDist(end,:) = [];
        massDist = [massDist; [1000*bodyDensity*bodyLength*pi*0.25*(bodyOuterDiameter^2 - bodyInnerDiameter^2),100*bodyLength,100*(noseLength + bodyLength/2)]];
        massDist(massDist(:,3)>100*(noseLength+bodyLength),:) = [];
        massCG = massDist(:,1)/1000;
        distCG = massDist(:,3)/100;
        rocketCG = sum(distCG.*massCG)/sum(massCG);
    end
    %</Stage separation>
    
    %<Interpolate for thrust, mass, and mass flow>
    if t<=max(timeStamp)
        thrustForce = interp1(timeStamp,Thrust,t,'linear');
        mass = interp1(timeStamp,Mass,t,'linear');
        massDot = interp1(timeStamp,MassDot,t,'linear');
    else
        thrustForce = 0;
        mass = Mass(end);
        massDot = 0;
    end
    %</Interpolate for thrust, mass, and mass flow>

    %<Numerical integral for velocity and height>
    dv = dt*(thrustForce - (0.5*rho*(0.007*0.145*4 + pi*(0.0792/2)^2)*Cd*v(end)^2) - mass*9.807)/mass;
    a= [a dv/dt];
    v = [v v(end)+dv];
    Mach = [Mach v(end)/sqrt(gamma*R*T)];

    h = [h h(end)+0.5*dt*(v(end)+v(end-1))];
    %</Numerical integral for velocity and height>

    %<Compute stability of rocket based on Mach and geometry>
    if Mach(end)<1
        [~,~,~,~,stability(end+1)] = subsonicRocket(bodyOuterDiameter,noseLength,bodyLength,numFinSets,Mach(end),0.01,finDims,rocketCG);
    else
        [~,~,~,~,stability(end+1)] = supersonicRocket(bodyOuterDiameter,noseLength,bodyLength,numFinSets,finResolution,Mach(end),0.01,finDims,rocketCG);
    end
    %</Compute stability of rocket based on Mach and geometry>
    
    %<Update drag coefficient based on Mach>
    if Mach(end)<0.2
        Cd = DragCoefficient(bodyDims,finDims,massDot,0.2,gamma,R,T,rho);
    else
        Cd = DragCoefficient(bodyDims,finDims,massDot,Mach(end),gamma,R,T,rho);
    end
    %</Update drag coefficient based on Mach>
    
    %<Update temp. and density based on ISA model>
    T = T0 - B*h(end);
    rho = rho0*(T/T0)^(g/(R*B) - 1);
    %</Update temp. and density based on ISA model>
    
    t = t+dt;
end

figure(1);
plot(dt:dt:t,h/0.3048,'LineWidth',3); % Lower bound for x may need to be changed to 0 or dt
xlabel("Time (s)");
ylabel("Height (ft)");
xlim([0 t]);
ylim([0 inf]);
title("Height after T");

figure(2);
xlabel("Time (s)");
yyaxis left
plot(dt:dt:t,v/0.3048,'LineWidth',3); % Lower bound for x may need to be changed to 0 or dt
ylabel("Velocity (ft/s)");
ylim([0 inf]);
yyaxis right
plot(dt:dt:t,Mach,'LineWidth',3); % Lower bound for x may need to be changed to 0 or dt
ylabel("Mach Number");
ylim([0 inf]);
xlim([0 t]);
title("Velocity/Mach after T");

figure(3);
plot(dt:dt:t,a/0.3048,'LineWidth',3); % Lower bound for x may need to be changed to 0 or dt
xlabel("Time (s)");
ylabel("Acceleration (ft/s^2)");
xlim([0 t]);
ylim([-inf inf]);
title("Acceleration after T");

figure(4);
plot(dt:dt:t,stability,'LineWidth',3); % Lower bound for x may need to be changed to 0 or dt
xlabel("Time (s)");
ylabel("Stability (cal)");
xlim([0 t]);
ylim([-inf inf]);
title("Stability after T");