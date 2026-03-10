clear; close all; clc; warning('off'); 
addpath('..'); addpath('../../Preliminary Files/Aerodynamics/Barrowman Functions');

simMaster = readcell("ZephyrMK3_SimMasterFileFake.xlsx");
setup = cell2mat(simMaster(2:11,2));
bodyDims = cell2mat(simMaster(2:8,4));
numStages = bodyDims(7);
bodyDims = cell2mat(simMaster(2:8+numStages,4));
finDims = cell2mat(simMaster(3:11,6:5+numStages));
massDist = cell2mat(simMaster(3:20,numStages+(7:9)));

gamma = setup(1);
R = setup(2);
p0 = setup(3)*1e5;
T0 = setup(4) + 273.15;
rho0 = p0/(R*T0);

bodyOuterDiameter = bodyDims(1)/100;
bodyInnerDiameter = bodyDims(2)/100;
noseConeLength = bodyDims(3)/100;
bodyDensity = bodyDims(4)*1000;
finDensity = bodyDims(5)*1000;
finResolution = bodyDims(6);
numStages = bodyDims(7);
bodyLength = sum(bodyDims(8:7+numStages))/100;

leadSweep = atand(finDims(9,:)./finDims(4,:));
trailSweep = atand(tand(leadSweep) - (finDims(2,:)-finDims(3,:))./finDims(4,:));
leadWedge = 2*atand(0.5*finDims(5,:)./finDims(6,:));
trailWedge = 2*atand(0.5*finDims(5,:)./finDims(7,:));
finDims = [finDims(1:8,:);leadSweep;leadSweep;trailSweep;trailSweep;leadWedge;trailWedge];

massDist = [massDist; [1e-3*bodyDensity*bodyLength*pi*0.25*(bodyOuterDiameter^2 - bodyInnerDiameter^2),100*bodyLength,100*(noseConeLength + bodyLength/2)];...
                      [1e-3*finDensity*0.5*finDims(4,:).*(finDims(2,:)+finDims(3,:)).*finDims(5,:);finDims(2,:);finDims(2,:)/2]'];
massDist(massDist(:,3)>100*(noseConeLength+bodyLength),:) = [];
massCG = massDist(:,1)/1000;
distCG = massDist(:,3)/100;
rocketCG = sum(distCG.*massCG)/sum(massCG);

motorThrustMass = uigetfile("*.xlsx");
motorThrustMass = readcell(motorThrustMass);
motort = cell2mat(motorThrustMass(2:end,1));
motorT = cell2mat(motorThrustMass(2:end,2));
motorM = cell2mat(motorThrustMass(2:end,3));

T = T0; p = p0; rho = rho0;
h = 0; v = 0; a = 0; Mach = 0; 
[~,~,~,~,stability] = subsonicRocket(bodyOuterDiameter,noseConeLength,bodyLength,numStages,Mach,0.01,finDims,rocketCG);
B = 0.0065; g = 9.807; Cd = 0; dt = 0.001; t = 0;
endCond = 'apogee'; endCondMet = 0;

while t<=max(motort)||~endCondMet
    t = t+dt
    %<Stage separation>
    if t == round(max(motort)+0.5,3)
        bodyDims(end) = [];
        finDims(:,end) = [];
        numStages = numStages - 1;
        bodyLength = sum(bodyDims(7:6+numStages))/100;
        massDist((end-numStages-1):end,:) = [];
        massDist = [massDist; [1e-3*bodyDensity*bodyLength*pi*0.25*(bodyOuterDiameter^2 - bodyInnerDiameter^2),100*bodyLength,100*(noseConeLength + bodyLength/2)];...
                              [1e-3*finDensity*0.5*finDims(4,:).*(finDims(2,:)+finDims(3,:)).*finDims(5,:);finDims(2,:);finDims(8,:)+finDims(2,:)/2]'];
        massDist(massDist(:,3)>100*(noseLength+bodyLength),:) = [];
        massCG = massDist(:,1)/1000;
        distCG = massDist(:,3)/100;
        rocketCG = sum(distCG.*massCG)/sum(massCG);
    end
    %</Stage separation>
    
    %<Interpolate for thrust, mass, and mass flow>
    if t<=max(motort)
        thrustForce = interp1(motort,motorT,t,'pchip');
        mass = interp1(motort,motorM,t,'pchip') + sum(massCG);
        massDot = (mass - interp1(motort,motorM,t-dt,'pchip'))/dt;
    else
        thrustForce = 0;
        mass = sum(massCG);
        massDot = 0;
    end
    %</Interpolate for thrust, mass, and mass flow>

    %<Numerical integral for velocity and height>
    dv = (dt/mass)*(thrustForce - mass*g - (0.5*rho*(v(end)^2)*Cd(end)*(4*max((finDims(4,:)/100).*(finDims(5,:)/100)) + 0.25*pi*bodyOuterDiameter^2)));
    a = [a dv/dt];
    v = [v v(end)+dv]; Mach = [Mach v(end)/sqrt(gamma*R*T)];
    h = [h h(end)+0.5*dt*(v(end)+v(end-1))];
    %</Numerical integral for velocity and height>

    %<Compute stability of rocket based on Mach and geometry>
    if Mach(end)<1
        [~,~,~,~,stability(end+1)] = subsonicRocket(bodyOuterDiameter,noseConeLength,bodyLength,numStages,Mach(end),0.01,finDims,rocketCG);
    else
        [~,~,~,~,stability(end+1)] = supersonicRocket(bodyOuterDiameter,noseConeLength,bodyLength,numStages,finResolution,Mach(end),0.01,finDims,rocketCG);
    end
    %</Compute stability of rocket based on Mach and geometry>
    
    %<Update drag coefficient based on Mach>
    if Mach(end)<0.2
        Cd = [Cd,DragCoefficient(bodyDims,finDims,massDot,0.2,gamma,p)];
    else
        Cd = [Cd,DragCoefficient(bodyDims,finDims,massDot,Mach(end),gamma,p)];
    end
    %</Update drag coefficient based on Mach>
    
    %<Update temp., pressure, and density based on ISA model>
    T = T0 - B*h(end);
    p = p0*(T/T0)^(g/(R*B));
    rho = p/(R*T);
    %</Update temp., pressure, and density based on ISA model>
    
    switch(endCond)
        case 'burnout'
            endCondMet = (t >= max(motort));
            %disp(100*t/max(motort));
        case 'separation'
            endCondMet = (t >= max(motort)+0.5);
            %disp(100*t/(max(motort)+0.5));
        case 'apogee'
            endCondMet = v(end)<=0;
            % if t<max(motort) disp(100*t/(max(motort)+0.5));
            % else disp(max(0,100*t/(t - v(end)/a(end))));
            % end
        case 'groundHit'
            endCondMet = h(end)<=0;
            % if t<max(motort) disp(100*t/(max(motort)+0.5));
            % else disp(max(0,real(100*t/(-v(end)/a(end) + sqrt((v(end)/a(end))^2 - 2*h(end)/a(end))))));
            % end
        otherwise
            error('Invalid simulaiton end condition');
    end
end

figure(1);
title("Height after T");
plot(0:dt:t,h/0.3048,'LineWidth',3); % Lower bound for x may need to be changed to 0 or dt
xlabel("Time (s)"); ylabel("Height (ft)");
xlim([0 t]); ylim([0 inf]);
saveas(figure(1),'Trajectory_Height3','png');

figure(2);
title("Velocity/Mach after T");
xlabel("Time (s)"); xlim([0 t]);
yyaxis left;
plot(0:dt:t,v/0.3048,'LineWidth',3); % Lower bound for x may need to be changed to 0 or dt
ylabel("Velocity (ft/s)"); ylim([0 inf]);
yyaxis right;
plot(0:dt:t,Mach,'LineWidth',3); % Lower bound for x may need to be changed to 0 or dt
ylabel("Mach Number"); ylim([0 inf]);
saveas(figure(2),'Trajectory_Velocity3','png');

figure(3);
title("Acceleration after T");
plot(0:dt:t,a/0.3048,'LineWidth',3); % Lower bound for x may need to be changed to 0 or dt
xlabel("Time (s)"); ylabel("Acceleration (ft/s^2)");
xlim([0 t]); ylim([-inf inf]);
saveas(figure(3),'Trajectory_Acceleration3','png');

figure(4);
title("Stability after T");
plot(0:dt:t,stability,'LineWidth',3); % Lower bound for x may need to be changed to 0 or dt
xlabel("Time (s)"); ylabel("Stability (cal)");
xlim([0 t]); ylim([1.5 8]);
saveas(figure(4),'Trajectory_Stability3','png');