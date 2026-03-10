clear; clc; close all; addpath('..');

gamma = 1.4; 
p_inf = 101325; %Pascals at MSL

simMaster = readcell("ZephyrMK2_SimMasterFile.xlsx");
setup = cell2mat(simMaster(2:9,2));
bodyDims = cell2mat(simMaster(2:7,4));
bodyDims = cell2mat(simMaster(2:7+bodyDims(6),4));
finDims = cell2mat(simMaster(2:15,6:6+bodyDims(6)-1));

machMin = setup(1);
machMax = setup(2);
machDelta = setup(3);
attackAngMin = setup(4)*pi/180;
attackAngMax = setup(5)*pi/180;
attackAngDelta = setup(6)*pi/180;
numel(attackAngMin:attackAngDelta:attackAngMax)

%Alocation
Cp1 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
p1 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
N1 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
A1 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
Cn1 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
Ca1 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);

Cp2 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
p2 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
N2 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
A2 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
Cn2 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
Ca2 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);

Cp3 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
p3 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
N3 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
A3 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
Cn3 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
Ca3 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);

Cp4 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
p4 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
N4 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
A4 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
Cn4 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
Ca4 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);

Cp5 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
p5 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
N5 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
A5 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
Cn5 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
Ca5 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);

Cp6 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
p6 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
N6 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
A6 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
Cn6 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
Ca6 = zeros(numel((machMin+machDelta):machDelta:machMax),numel(attackAngMin:attackAngDelta:attackAngMax),2);
%%
for n = 1:bodyDims(6)

    TaperLength = finDims(6,n)/100;
    RootChord = finDims(2,n)/100;
    MidLength = RootChord - 2*TaperLength;
    t = finDims(5,n)/100;
    EdgeLength = sqrt((TaperLength^2)+((t/2)^2)); 
    Lphi = (finDims(13,n)/2)*pi/180;
    Tphi = (finDims(14,n)/2)*pi/180;

    j = 1;
    for aoa = attackAngMin:attackAngDelta:attackAngMax

       i = 1;
         for Mach = (machMin+machDelta):machDelta:machMax

            q_inf = (gamma/2)*p_inf*(Mach^2);

            Cp1(i,j,n) = 2*(Lphi-aoa)/sqrt((Mach^2)-1); %Top leading edge
            p1(i,j,n) = p_inf+(q_inf*Cp1(i,j,n));
            N1(i,j,n) = -p1(i,j,n)*EdgeLength*cos(Lphi);
            A1(i,j,n) = p1(i,j,n)*EdgeLength*sin(Lphi);
            Cn1(i,j,n) = N1(i,j,n)/(q_inf*RootChord);
            Ca1(i,j,n) = A1(i,j,n)/(q_inf*RootChord);

            Cp2(i,j,n) = 2*(-aoa)/sqrt((Mach^2)-1); %Top middle section
            p2(i,j,n) = p_inf+(q_inf*Cp2(i,j,n));
            N2(i,j,n) = -p2(i,j,n)*MidLength;
            A2(i,j,n) = 0;
            Cn2(i,j,n) = N2(i,j,n)/(q_inf*RootChord);
            Ca2(i,j,n) = A2(i,j,n)/(q_inf*RootChord);

            Cp3(i,j,n) = 2*(-Tphi-aoa)/sqrt((Mach^2)-1); %Top trailing edge
            p3(i,j,n) = p_inf+(q_inf*Cp3(i,j,n));
            N3(i,j,n) = -p3(i,j,n)*EdgeLength*cos(Tphi);
            A3(i,j,n) = -p3(i,j,n)*EdgeLength*sin(Tphi);
            Cn3(i,j,n) = N3(i,j,n)/(q_inf*RootChord);
            Ca3(i,j,n) = A3(i,j,n)/(q_inf*RootChord);
    
            Cp4(i,j,n) = 2*(Lphi+aoa)/sqrt((Mach^2)-1); %Bottom leading edge
            p4(i,j,n) = p_inf+(q_inf*Cp4(i,j,n));
            N4(i,j,n) = p4(i,j,n)*EdgeLength*cos(Lphi);
            A4(i,j,n) = p4(i,j,n)*EdgeLength*sin(Lphi);
            Cn4(i,j,n) = N4(i,j,n)/(q_inf*RootChord);
            Ca4(i,j,n) = A4(i,j,n)/(q_inf*RootChord);
    
            Cp5(i,j,n) = 2*(aoa)/sqrt((Mach^2)-1); %Bottom middle section
            p5(i,j,n) = p_inf+(q_inf*Cp5(i,j,n));
            N5(i,j,n) = p5(i,j,n)*MidLength;
            A5(i,j,n) = 0;
            Cn5(i,j,n) = N5(i,j,n)/(q_inf*RootChord);
            Ca5(i,j,n) = A5(i,j,n)/(q_inf*RootChord);
            
            Cp6(i,j,n) = 2*(-Tphi+aoa)/sqrt((Mach^2)-1); %Bottom trailing edge
            p6(i,j,n) = p_inf+(q_inf*Cp6(i,j,n));
            N6(i,j,n) = p6(i,j,n)*EdgeLength*cos(Tphi);
            A6(i,j,n) = -p6(i,j,n)*EdgeLength*sin(Tphi);
            Cn6(i,j,n) = N6(i,j,n)/(q_inf*RootChord);
            Ca6(i,j,n) = A6(i,j,n)/(q_inf*RootChord);
    

            i = i + 1;
         end
         j = j + 1;

    end
end

Cn_sup = Cn1 + Cn2 + Cn3 + Cn4 + Cn5 + Cn6;
Ca_sup = Ca1 + Ca2 + Ca3 + Ca4 + Ca5 + Ca6;

save('Supersonic_Axial','Ca_sup');
save('Supersonic_Normal','Cn_sup');