function [CLa,centerPresX,liftForce,dragForce] = subsonicFin(bodyRadius,finDims,Mach,attackAng)

numFins = finDims(1); rootChord = finDims(2)/100; 
tipChord = finDims(3)/100; semispan = finDims(4)/100; finThickness = finDims(5)/100;
distNoseTip = finDims(8)/100; leadSweep = finDims(9)*pi/180;  

t = (semispan + bodyRadius)/(bodyRadius);
KTB = (2/(pi*(1 - t^-1)^2))*((1 + t^-4)*(0.5*atan(0.5*(t - t^-1))) - (t^-2)*(t - t^-1 + 2*atan(1/t)));
KBT = ((1 + 1/t)^2) - KTB;

midSweep = atan(tan(leadSweep) - 0.5*(rootChord - tipChord)/semispan);

finArea = 0.5*(rootChord + tipChord) * semispan;
AspectRatio = (semispan)^2 / finArea;

refArea = pi*(bodyRadius)^2;
CLaT = (numFins* pi * (AspectRatio * finArea / refArea)) / (2 + sqrt(4 + ( (sqrt(abs(1 - Mach^2))) * AspectRatio / cos(midSweep) )^2 ));

xT = semispan * tan(leadSweep);
xBarT = distNoseTip + (xT/3) * ((rootChord + 2 * tipChord) / (rootChord + tipChord)) + (1/6) * (rootChord + tipChord - (rootChord * tipChord) / (rootChord + tipChord)) ;

CLa(1) = CLaT * KTB;
CLa(2) = CLaT * KBT;

centerPresX(1) = xBarT;
centerPresX(2) = distNoseTip + rootChord/4 + ((sqrt(semispan^2 - bodyRadius^2)*acosh(semispan/bodyRadius) - semispan + pi*bodyRadius/2)/((bodyRadius/sqrt(semispan^2 - bodyRadius^2))*acosh(semispan/bodyRadius) + semispan/bodyRadius - pi/2) - bodyRadius)*(tan(leadSweep) - 0.25*(rootChord-tipChord)/semispan);

liftForce = (0.5*1.225*(Mach*340.3)^2)*CLa*attackAng*refArea;

cMAC = (2/3)*(rootChord + tipChord - (rootChord*tipChord)/(rootChord + tipChord));
Re = 1.225*(Mach*340.3)*cMAC/(1.789*10^-5);

wettedArea = 2*numFins*finArea;

Cfc = (1 - 0.12*Mach^2)*((3.46*log(Re) - 5.6)^-2);
CDf = Cfc*wettedArea/refArea;

if Mach<0.9, dCD = (1-Mach^2)^(-0.417) - 1; else, dCD = 1 - 1.5*(Mach-0.9); end

CDL = 2*numFins*(semispan*0.01/refArea)*(cos(leadSweep)^2)*dCD;

Cfc_1 = 0.88*(3.46*log(Re) - 5.6)^-2;
K = cos(midSweep)^2 + ((0.223 + 4.02*Cfc_1)^2)/(((rootChord/finThickness)*Cfc_1)^(2/3));
Pf = (K - (Mach*cos(midSweep))^2)^(-1/2);
CfB = 2*((3.46*log(Re) - 5.6)^-2)*(rootChord/finThickness);

%%% CDB = 0.135*numFins*(finArea/refArea)*Pf*CfB^(-1/3);

Cf_laminar = 1.328/(sqrt(Re));
Cf_turbulent = (3.46*log(Re) - 5.6)^(-2) - 1700/Re;

if Re>=500000 % turbulent
    Cfc = Cf_turbulent * (1 - 1.2*Mach^2);
    Cfc_sonic = Cf_turbulent * (1 - 1.2*1^2);
else % laminar
    Cfc = Cf_laminar * (1 - 0.9*Mach^2);
    Cfc_sonic = Cf_laminar * (1 - 0.9*1^2);
end

% K_cdb = cos(Gamma_c)^2 + [(0.223 + 4.02*Cfc_sonic*(finThickness/finThickness)^2)*(Cfc_sonic*(rootChord/finThickness))^(-1/3)]^2
K_cdb = cos(midSweep)^2 + ((0.223 + 4.02*Cfc_sonic*(finThickness/finThickness)^2)*(Cfc_sonic*(rootChord/finThickness))^(-1/3))^2;

p = (K_cdb - (Mach *cos(midSweep))^2)^(-1/2);

% CDB = 0.135*N*(A_bf/A_r)*((2*C_fc*(C_r/h_r))^(-1/3))*p
CDB = 0.135*numFins*(finArea/refArea)*((2*Cfc*(rootChord/finThickness))^(-1/3))*p;

% No CDP for the fins

xi = AspectRatio*(finThickness/rootChord)^(1/3);
CDT_1 = 1.15*((finThickness/rootChord)^(5/3))*(1.61 + xi - sqrt((xi - 1.43)^2 + 0.578));
K = cos(midSweep)^2 + (((CDT_1/Cfc_1)*(refArea/finArea) - 4*(finThickness/rootChord)*cos(midSweep))/(120*((finThickness/rootChord)^4)*cos(midSweep)^2))^(2/3);
Pf = (K - (Mach*cos(midSweep))^2)^(-1/2);
CDT = 4*numFins*Cfc*(finArea/refArea)*((finThickness/rootChord)*cos(midSweep) + 30*((finThickness/rootChord)^4)*(cos(midSweep)^2)*Pf^3);

CD = CDf + CDL + CDB + CDT;
dragForce = real(CD*(0.5*1.225*(Mach*340.3)^2)*finArea);