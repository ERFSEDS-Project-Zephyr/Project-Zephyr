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

liftForce = (0.5*1.225*(Mach*340.3)^2)*numFins*CLa*attackAng;

cMAC = (2/3)*(rootChord + tipChord - (rootChord*tipChord)/(rootChord + tipChord));
Re = 1.225*(Mach*340.3)*cMAC/(1.789*10^-5);
Cfc = (1.328/sqrt(Re))*(1 - 0.12*Mach^2);
CDF = Cfc*2*numFins*finArea/refArea;
chi = AspectRatio*(finThickness/rootChord)^(1/3);
cdtt = 1.15*((finThickness/rootChord)^(5/3))*(1.61 + chi - sqrt((chi - 1.43)^2 + 0.578));
K = cos(midSweep)^2 + ((0.223 + 4.02*Cfc)^2)/((Cfc*rootChord/finThickness)^(2/3));
CDB = 0.135*numFins/(((2*Cfc*rootChord/finThickness)^(1/3))*sqrt(K - (Mach*cos(midSweep))^2));
K = cos(midSweep)^2 + (((cdtt/Cfc)*(refArea/finArea) - 4*(finThickness/rootChord)*cos(midSweep))/(120*((finThickness/rootChord)^4)*cos(midSweep)^2))^(2/3);
CDT = 4*numFins*Cfc*((finThickness/rootChord)*cos(midSweep) + 30*((finThickness/rootChord)^4)*(cos(midSweep)^2)/sqrt((K - (Mach*cos(midSweep))^2)^3));

CD = CDF + CDB + CDT;
dragForce = CD*(0.5*1.225*(Mach*340.3)^2)*numFins*finArea/refArea;