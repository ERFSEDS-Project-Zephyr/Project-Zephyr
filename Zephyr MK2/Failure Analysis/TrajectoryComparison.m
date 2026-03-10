clear;
clc;
close all;
addpath('..');
addpath('../Aerodynamics');
addpath('../Aerodynamics/Barrowman Functions');

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
Cd = 0;
dt = 0.001;
t = dt;

g = 9.807;
B = 0.0065;
gamma = 1.4;
R = 8.314/0.0288;
%</Establish initial conditions and universal constants>

numFinSets = 2;

while t<=max(timeStamp)||h(end)>0
    %<Stage separation>
    if any(t == round(timeStamp(MassDot < -5),3))
        numStages = numStages - 1;
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
    a = [a (thrustForce - (0.5*rho*(0.007*0.145*4 + pi*(bodyOuterDiameter/2)^2)*Cd*v(end)^2) - mass*9.807)/mass];
    v = [v v(end)+0.5*dt*(a(end)+a(end-1))];
    Mach = [Mach v(end)/sqrt(gamma*R*T)];

    h = [h h(end)+0.5*dt*(v(end)+v(end-1))];
    %</Numerical integral for velocity and height>
    
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

dart = readcell('ZephyrMk2Flight1_Dart.csv');
booster = readcell('Zephyrmk2_flight1_booster.csv');

dart_time = cell2mat(dart(2:end,5)); T = find(dart_time==0,1,'last');
dart_time = dart_time(T:end);
dart_height = cell2mat(dart((T+1):end,11));
dart_accelx = cell2mat(dart((T+1):end,17));
dart_accely = cell2mat(dart((T+1):end,18));
dart_accelz = cell2mat(dart((T+1):end,19));
dart_accel = sign(dart_accelx).*sqrt(dart_accelx.^2 + dart_accely.^2 + dart_accelz.^2);
dart_veloc = 0;
for i = 2:length(dart_time)
    dart_veloc(i,:) = dart_veloc(i-1) + 0.5*(dart_accel(i) + dart_accel(i-1))*(dart_time(i) - dart_time(i-1));
end

booster_time = cell2mat(booster(2:end,5)); T = find(booster_time==0,1,'last');
booster_time = booster_time(T:end);
booster_height = cell2mat(booster((T+1):end,11));
booster_accelx = cell2mat(booster((T+1):end,17));
booster_accely = cell2mat(booster((T+1):end,18));
booster_accelz = cell2mat(booster((T+1):end,19));
booster_accel = sign(booster_accelx).*sqrt(booster_accelx.^2 + booster_accely.^2 + booster_accelz.^2);
booster_veloc = 0;
for i = 2:length(booster_time)
    booster_veloc(i,:) = booster_veloc(i-1) + 0.5*(booster_accel(i) + booster_accel(i-1))*(booster_time(i) - booster_time(i-1));
end

h_dart = interp1(0:dt:t,h,dart_time,'linear','extrap');
h_booster = interp1(0:dt:t,h,booster_time,'linear','extrap');
v_dart = interp1(0:dt:t,v,dart_time,'linear','extrap');
v_booster = interp1(0:dt:t,v,booster_time,'linear','extrap');
a_dart = interp1(0:dt:t,a,dart_time,'linear','extrap');
a_booster = interp1(0:dt:t,a,booster_time,'linear','extrap');

figure(1); hold on;
plot(dart_time,(dart_height-h_dart)./dart_height,'LineWidth',3);
plot(booster_time,(booster_height-h_booster)./booster_height,'LineWidth',3);
xlabel("Time (s)");
ylabel("Height (ft)");
xlim([0 t]); ylim([-1 1]);
title("Height after T");
legend('Dart','Booster','Orientation','horizontal','Location','southoutside');

figure(2); hold on;
plot(dart_time,(dart_veloc-v_dart)./dart_veloc,'LineWidth',3);
plot(booster_time,(booster_veloc-v_booster)./booster_veloc,'LineWidth',3);
xlabel("Time (s)");
ylabel("Velocity (ft/s)");
xlim([0 t]); ylim([-1 1]);
title("Velocity after T");
legend('Dart','Booster','Orientation','horizontal','Location','southoutside');

figure(3); hold on;
plot(dart_time,(dart_accel-a_dart)./dart_accel,'LineWidth',3);
plot(booster_time,(booster_accel-a_booster)./booster_accel,'LineWidth',3);
xlabel("Time (s)");
ylabel("Acceleration (ft/s^2)");
xlim([0 t]); ylim([-1 1]);
title("Acceleration after T");
legend('Dart','Booster','Orientation','horizontal','Location','southoutside');