function [CLa,centerPresX,liftForce,dragForce] = subsonicNose(bodyOuterDiameter,noseLength,bodyLength,Mach,attackAng)
bodyRadius = bodyOuterDiameter/2;
ogiveRadius = (noseLength^2 + bodyRadius^2)/(2*bodyRadius);
filledVolume = pi*((bodyRadius-ogiveRadius)*asin(noseLength/ogiveRadius)*ogiveRadius^2 + (bodyRadius^2 - 2*ogiveRadius*bodyRadius + 2*ogiveRadius^2)*noseLength + (bodyRadius-ogiveRadius)*noseLength*sqrt(ogiveRadius^2 - noseLength^2) - (noseLength^3)/3);
CLa = 2; % Eq. 3-66
centerPresX = noseLength - filledVolume/(pi*bodyRadius^2); % Eq. 3-89

fN = noseLength/bodyOuterDiameter; fB = bodyLength/bodyOuterDiameter; F = ogiveRadius/noseLength;
q = 0.5*1.225*(Mach*340.3)^2;
Anormal = (noseLength^2)*((F^2)*asin(1/F) - sqrt(F^2 - 1));
liftForce = (attackAng*CLa)*q*Anormal;

areaWetted = pi*(noseLength^2)*((F^2)*asin(1/F) + 3*sqrt((F^2) - 1));
Re = 1.225*(Mach*340.3)*noseLength/(1.789*10^-5);
Re_1 = 1.225*340.3*noseLength/(1.789*10^-5);

if Re < 500000 && Mach<1
    Cfc = (1.328/sqrt(Re))*(1 - 0.9*Mach^2);
    Cfc_1 = (1.328/sqrt(Re_1))*(0.1);
elseif Re < 500000 && Mach>1
    Cfc = ((3.46*log(Re) - 5.6)^-2 - 1700/Re)*(1 - 1.2*Mach^2);
    Cfc_1 = ((3.46*log(Re_1) - 5.6)^-2 - 1700/Re_1)*(-0.2);
elseif Re >= 500000 && Mach<1
    Cfc = (1.328/sqrt(Re))*(1 + 0.15*Mach^2)^(-1/4);
    Cfc_1 = (1.328/sqrt(Re_1))*(1.15^(-1/4));
else
    Cfc = ((3.46*log(Re) - 5.6)^-2 - 1700/Re)*(1 + 0.18*Mach^2)^-1;
    Cfc_1 = ((3.46*log(Re_1) - 5.6)^-2 - 1700/Re_1)*(1.18^-1);
end
CDf = (1 + 1/(2*fB))*Cfc*areaWetted;

K = 1 + ((6.82*Cfc_1*areaWetted*(fN+1)^2.22)/(fB^3))^(5/3);
CDP = 6*Cfc*areaWetted/((fB^3)*(K - Mach^2)^0.6);

CD = CDf + CDP;
if Mach==0
    dragForce = 0;
else
    dragForce = real((0.5*1.225*(Mach*340.3)^2)*CD*(pi*bodyRadius^2));
end