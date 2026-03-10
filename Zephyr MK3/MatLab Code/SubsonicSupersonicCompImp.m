%Clear workspace
clc
close all
clear
addpath('..\Data Files\')

%Read data from sim file
simMaster = readcell("ZephyrMK3_SimMasterFile.xlsx");
setup = cell2mat(simMaster(3:15,4));
bodyDims = cell2mat(simMaster(3:7,8));
% bodyDims = cell2mat(simMaster(2:8+bodyDims(7),4));
finDims = cell2mat(simMaster(3:12,12:13));

%Values unique to each fin set
chordRoot = finDims(2,:); % [in]
chordTip = finDims(3,:); % [in]
s = finDims(4,:); % [in]
t = finDims(6,:); % [in]

%%
%Atmospheric Values
gamma=setup(1);
GasConstant = setup(2); % [ft*lbf/(slug*Rankine)]
p_inf=setup(4); % [psi]
temperature_inf = setup(5)+ 459.67; % [F] to [Rankine]
mu_inf = 3.62*10^-7; % [slug/(ft*s)] Viscosity of air at SSL by NASA
a = sqrt(gamma*GasConstant*temperature_inf); %Speed of sound at SSL, [ft/s]
rho_inf = (p_inf*144)/(GasConstant*temperature_inf); % [slug/ft^3] of air

%Test Parameters
minMach=setup(6);
maxMach=setup(7);
deltaMach=setup(8);
minAngleAttack=setup(9);
maxAngleAttack=setup(10);
DeltaAngleAttack=setup(11);


%Initializing Tables
MachVect = minMach:deltaMach:maxMach;
alphaVect = minAngleAttack:DeltaAngleAttack:maxAngleAttack;
NormForceTable = zeros(length(MachVect),length(alphaVect)); 
AxialForceTable = zeros(length(MachVect),length(alphaVect)); 
k = 1;

for finset = 1:2

    sTemp = s(finset); % [in]

    j = 1;

    for AngleAttack = deg2rad((minAngleAttack:DeltaAngleAttack:maxAngleAttack))
        
        i = 1;
        
        %Subsonic
        for Mach = minMach+deltaMach:deltaMach:1-deltaMach

            c_D=(2*pi*AngleAttack^2)/sqrt(1-(Mach^2)); 
            c_L=(2*pi*AngleAttack)/sqrt(1-(Mach^2)); 
            p_stream=p_inf*(1+(gamma-1/2)*Mach^2)^(-gamma/(gamma-1)); % [psi]
            q_inf=(gamma/2)*p_stream*(Mach)^2; % [psi]

            D(i,j,finset)=c_D*q_inf*sTemp; % [lb/in]
            L(i,j,finset)=c_L*q_inf*sTemp; % [lb/in]
            NormForce = L*cos(AngleAttack)+D*sin(AngleAttack); % [lb/in]
            Axial = D*cos(AngleAttack)+L*sin(AngleAttack); % [lb/in]
        
                
            Ca(i,j,finset) = (c_D.*cos(AngleAttack) + c_L.*sin(AngleAttack));
            Cn(i,j,finset) = (-c_D.*sin(AngleAttack) + c_L.*cos(AngleAttack));
                
            i=i+1;
        end

        %Supersonic
        for Mach = 1+deltaMach:deltaMach:maxMach
        
            c_D=(4*AngleAttack)/sqrt(Mach^2-1);
            c_L=(4*AngleAttack)/sqrt(Mach^2-1);
            p_stream=p_inf*(1+(gamma-1/2)*Mach^2)^(-gamma/(gamma-1)); % [psi]
            q_inf=(gamma/2)*p_stream*(Mach)^2; % [psi]

            D(i,j,finset)=c_D*q_inf*sTemp; % [lb/in]
            L(i,j,finset)=c_L*q_inf*sTemp; % [lb/in]
            NormForce = L*cos(AngleAttack)+D*sin(AngleAttack); % [lb/in]
            Axial = D*cos(AngleAttack)+L*sin(AngleAttack); % [lb/in]
        
                
            Ca(i,j,finset) = c_D.*cos(AngleAttack) + c_L.*sin(AngleAttack);
            Cn(i,j,finset) = -c_D.*sin(AngleAttack) + c_L.*cos(AngleAttack);

            i = i + 1;

        end

        j = j + 1;
        k = k + 1;

    end
end


%% CD calculations for other parts @ AngleAttack  = 0

ConeHeight = bodyDims(3); % [in]
TubeRadius = bodyDims(1)/2; % radius is half diameter [in]
WedgeAngle = atan(TubeRadius/ConeHeight); % assume triangle cone for TAT [rad]
TubeLength = finDims(9,:)+finDims(2,:)-ConeHeight; % [in]

Cnose = sqrt(ConeHeight^2+TubeRadius^2); % length of nose assume triangle [in]

rho_SSL = rho_inf; % [slug/ft^3] of air


CDOindex=1; % each column is increasing mach, not dependent on angle
for Mach = 0+deltaMach:deltaMach:maxMach
    %Isentropic relations
    rho(CDOindex) = rho_inf*(1 + (gamma-1)/2*Mach^2)^(-1/(gamma-1)); % [slug/ft^3]
    temperature(CDOindex) = temperature_inf*(1 + (gamma-1)/2*Mach^2)^(-1); % [Rankine]
    % Sutherland for viscosity
    mu(CDOindex) = mu_inf*(temperature(CDOindex)/temperature_inf)^(3/2)*(temperature_inf+198.72)/(temperature(CDOindex)+198.72);
        % [slug/(ft*s)]
    for k = 1:2
        ReFins(k) = (rho_inf*Mach*a*((chordRoot(k) + chordTip(k))/2))/mu(CDOindex);
        CDOFins(k,CDOindex) = (0.074/(ReFins(k)^(1/5)))*8*s(k);
    end
    
    Vnose = Mach*a/cos(0.5*WedgeAngle); % [ft/s] idk if half is needed
    ReNose = (rho_inf*Mach*a*Cnose)/mu(CDOindex);
    ReBody = (rho_inf*Mach*a*(TubeLength(1)+TubeLength(2)))/mu(CDOindex);

    CDOBody(CDOindex) = (0.074/(ReBody^(1/5)))*2*pi*TubeRadius;
    CDONose(CDOindex) = (0.074/(ReNose^(1/5)))*pi*TubeRadius;
    CDOindex = CDOindex + 1;
end
% CD0 is like skinfriction drag so added to all entries
Ca(:,:,1) = Ca(:,:,1) + CDOFins(1,[1:99,101:end])'*ones(1,size(Ca,2)); %skipping mach 1
Ca(:,:,2) = Ca(:,:,2) + CDOFins(2,[1:99,101:end])'*ones(1,size(Ca,2));

%% 2D Graph of axial Coeff
% Mach index skiping 1
MachIndex = [(minMach+deltaMach):deltaMach:(1-deltaMach),(1+deltaMach):deltaMach:maxMach];
AngleIndex = (minAngleAttack:DeltaAngleAttack:maxAngleAttack);
%Finset 1
figure(Name="2D Axial Coeff Finset 1",NumberTitle="off");
hold on
plot(MachIndex,Ca(:,:,1))
xlabel("Mach [M]")
ylabel("Coeff Axial")
title("Coeff Axial vs Mach Finset 1")
hold off
%Finset 2 (redundant)
%{
figure(Name="2D Axial Coeff Finset 2",NumberTitle="off");
hold on
plot(MachIndex,Ca(:,:,2))
xlabel("Mach [M]")
ylabel("Coeff Axial")
title("Coeff Axial vs Mach Finset 2")
hold off
%}
%% 3D Surf Graph of Axial Coeff

[Y, X] = meshgrid(AngleIndex,MachIndex);
%Finset 1
figure(Name="3D Axial Coeff Finset 1",NumberTitle="off");
hold on
AxialGraph1 = surf(X,Y,Ca(:,:,1));
grid on
ylabel("Angle Of Attack")
xlabel("Mach Number")
zlabel("Coefficient of Axial Force")
title("Axial Coeff Finset 1")
view([-30, 30]);
hold off
%Finset 2 (redundant) 
%{
figure(Name="3D Axial Coeff Finset 2",NumberTitle="off");
hold on
AxialGraph2 = surf(X,Y,Ca(:,:,2));
grid on
ylabel("Angle Of Attack")
xlabel("Mach Number")
zlabel("Coefficient of Axial Force")
title("Axial Coeff Finset 2")
view([-30, 30]);
hold off
%}
%% 2D Graph of Lift and Drag
%Lift

%Finset 1
figure(Name="2D Lift Finset 1",NumberTitle="off");
hold on
plot(MachIndex,L(:,:,1))
xlabel("Mach [M]")
ylabel("Lift")
title("Lift Finset 1")
hold off
%Finset 2
figure(Name="2D Lift Finset 2",NumberTitle="off");
hold on
plot(MachIndex,L(:,:,2))
xlabel("Mach [M]")
ylabel("Lift")
title("Lift Finset 2")
hold off

%Drag

%Finset 1
figure(Name="2D Drag Finset 1",NumberTitle="off");
hold on
plot(MachIndex,D(:,:,1))
xlabel("Mach [M]")
ylabel("Drag")
title("Drag Finset 1")
hold off
%Finset 2
figure(Name="2D Drag Finset 2",NumberTitle="off");
hold on
plot(MachIndex,D(:,:,2))
xlabel("Mach [M]")
ylabel("Drag")
title("Drag Finset 2")
hold off


%% 3D Surf Graphs of Lift and Drag
%Lift

%Finset 1
figure(Name="3D Lift Finset 1",NumberTitle="off");
hold on
LiftGraph1 = surf(X,Y,L(:,:,1));
grid on
ylabel("Angle Of Attack")
xlabel("Mach Number")
zlabel("Lift")
title("Lift Finset 1")
view([-30, 30]);
hold off
%Finset 2
figure(Name="3D Lift Finset 2",NumberTitle="off");
hold on
LiftGraph2 = surf(X,Y,L(:,:,2));
grid on
ylabel("Angle Of Attack")
xlabel("Mach Number")
zlabel("Lift")
title("Lift Finset 2")
view([-30, 30]);
hold off

%Drag

%Finset 1
figure(Name="3D Drag Finset 1",NumberTitle="off");
hold on
DragGraph1 = surf(X,Y,D(:,:,1));
grid on
ylabel("Angle Of Attack")
xlabel("Mach Number")
zlabel("Drag")
title("Drag Finset 1")
view([-30, 30]);
hold off
%Finset 2
figure(Name="3D Drag Finset 2",NumberTitle="off");
hold on
DragGraph2 = surf(X,Y,D(:,:,2));
grid on
ylabel("Angle Of Attack")
xlabel("Mach Number")
zlabel("Drag")
title("Drag Finset 2")
view([-30, 30]);
hold off
