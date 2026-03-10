function [CLa,centerPresX,totalLift] = supersonicFin(bodyDiameter,finResolution,Mach,attackAng,finDims,k)
liftForce = 0;
pitchMoment = 0;

numFins = finDims(1); rootChord = finDims(2)/100; tipChord = finDims(3)/100; semispan = finDims(4)/100; leadTaper = finDims(5)/100; trailTaper = finDims(6)/100; distNoseTip = finDims(7)/100;
leadSweep = finDims(8)*pi/180; leadTaperSweep = finDims(9)*pi/180; trailTaperSweep = finDims(10)*pi/180; trailSweep = finDims(11)*pi/180; leadWedge = finDims(12)*pi/360; trailWedge = finDims(13)*pi/360;

refArea = pi*(bodyDiameter^2)/4;
dy = semispan/finResolution;
coeffScale = 2*dy/refArea;

g = 1.4;
B = sqrt(abs(1 - Mach^2));
K1 = 2/B; % Eq. 2a
K2 = ((g+1)*Mach^4 - 4*B^2)/(2*B^4); % Eq. 2b
K3 = ((g+1)*Mach^8 + (2*g^2 - 7*g - 5)*Mach^6 + 10*(g+1)*Mach^4 - 12*Mach^2 + 8)/(6*B^7); % Eq. 2b
Kstar = (g+1)*(Mach^4)*((5 - 3*g)*Mach^4 + 4*(g-3)*(Mach^2) + 8)/(48*B^7); % Eq. 3

incl(1) = leadWedge - attackAng;
incl(2) = -attackAng;
incl(3) = -trailWedge - attackAng;
incl(4) = leadWedge + attackAng;
incl(5) = attackAng;
incl(6) = -trailWedge + attackAng;

for I=1:6
    regionPresCoeff(I) = K1*incl(I) + K2*incl(I)^2 + K3*incl(I)^3;
    if I>=4
        regionPresCoeff(I) = regionPresCoeff(I) - Kstar*incl(4)^3;
    elseif incl(1) > 0
        regionPresCoeff(I) = regionPresCoeff(I) - Kstar*incl(1)^3;
    else
        continue;
    end
end

y = 0;
while y <= semispan
    sectionChord = y*(tan(trailSweep) - tan(leadSweep)) + rootChord;
    sectionLeadChord = y*(tan(leadTaperSweep) - tan(leadSweep)) + leadTaper;
    sectionTrailChord = y*(tan(trailSweep) - tan(trailTaperSweep)) + trailTaper;
    sectionLeadLength = sectionLeadChord/cos(leadWedge);
    sectionMidLength = sectionChord - sectionLeadChord - sectionTrailChord;
    sectionTrailLength = sectionTrailChord/cos(trailWedge);
    distLeadRoot = y*tan(leadSweep)*cos(attackAng);
    regionStreamLength(1) = sectionLeadLength*cos(incl(1));
    regionStreamLength(2) = sectionMidLength*cos(incl(2));
    regionStreamLength(3) = sectionTrailLength*cos(incl(3));
    regionStreamLength(4) = sectionLeadLength*cos(incl(4));
    regionStreamLength(5) = sectionMidLength*cos(incl(5));
    regionStreamLength(6) = sectionTrailLength*cos(incl(6));

    machConeRatio = [0,0,0,0,0,0];

    if k==1
        LW = (semispan-y)*(tan(leadSweep) + B);
        if LW >= sectionLeadLength
            machConeRatio(3) = (sectionChord-LW)/sectionTrailChord;
            machConeRatio(6) = machConeRatio(3);
        elseif LW >=sectionLeadChord
            machConeRatio(3) = 1;
            machConeRatio(6) = machConeRatio(3);
            machConeRatio(2) = (sectionChord-sectionTrailChord-LW)/sectionMidLength;
             machConeRatio(5) = machConeRatio(2);
        else
            machConeRatio(3) = 1;
            machConeRatio(6) = machConeRatio(3);
            machConeRatio(2) = 1;
            machConeRatio(5) = machConeRatio(2);
            machConeRatio(1) = (sectionLeadChord-LW)/sectionLeadChord;
            machConeRatio(4) = machConeRatio(1);
        end
    end

    for I=1:6
        if I<=3
            regionLiftForce(I) = -regionPresCoeff(I)*regionStreamLength(I)*(1 - 0.5*machConeRatio(I));
        else
            regionLiftForce(I) = regionPresCoeff(I)*regionStreamLength(I)*(1 - 0.5*machConeRatio(I));
        end
    end
    sectionEnd(1) = distLeadRoot;
    sectionEnd(2) = distLeadRoot + regionStreamLength(1);
    sectionEnd(3) = distLeadRoot + regionStreamLength(1) + regionStreamLength(2);
    sectionEnd(4) = distLeadRoot;
    sectionEnd(5) = distLeadRoot + regionStreamLength(4);
    sectionEnd(6) = distLeadRoot + regionStreamLength(4) + regionStreamLength(5);

    sectionPitchMoment = 0;
    for I=1:6
        if regionStreamLength(I)~=0
            regionMomentArm = 0.5*regionStreamLength(I)*(1-machConeRatio(I)+0.5*machConeRatio(I)^2 + sectionEnd(I)*(2-machConeRatio(I))/regionStreamLength(I))/(1-0.5*machConeRatio(I));
        else
            regionMomentArm = 0;
        end
        sectionPitchMoment = sectionPitchMoment - regionMomentArm*regionLiftForce(I);
    end
    sectionLiftForce = 0;
    for I=1:6
        sectionLiftForce = sectionLiftForce + regionLiftForce(I);
    end
    liftForce = liftForce + sectionLiftForce;
    pitchMoment = pitchMoment + sectionPitchMoment;
    y = y + dy;
end

bodyRadius = bodyDiameter/2;
t = (semispan + bodyRadius)/bodyRadius;
s = t - 1;
m = cot(leadSweep);
CLaT = (180/pi)*numFins*liftForce*coeffScale;

KTB = (2/(pi*(1 - t^-1)^2))*((1 + t^-4)*(0.5*atan(0.5*(t - t^-1))) - (t^-2)*(t - t^-1 + 2*atan(1/t)));
KBT = (1 + t^-2)^2 - KTB;
% if Mach <= sec(leadSweep)
%     KBT = ((1 - s^-2) - (2/pi)*((1 + s^4)*(0.5*atan(0.5*(s - s^-1)) + pi/4) - (s^-2)*(s - s^-1 + 2*atan(s^-1))))*(1 - s^-1)^-2;
% else
%     KBT = ((8*m)/(pi*sqrt((B*m)^2 - 1)*(1 + tipChord/rootChord)*(B*bodyDiameter/rootChord)*(t-2)*CLaT))*((B*m/(1 + B*m))*((((B*m + 1)*(B*bodyDiameter/rootChord) + B*m)/(B*m))^2)*acos((1 + (1+B*m)*B*bodyDiameter/rootChord)/(B*m + (1+B*m)*B*bodyDiameter/rootChord) +...
%                                                                                                              (sqrt((B*m)^2 - 1)/(B*m + 1))*(sqrt(1 + 2*B*bodyDiameter/rootChord) - 1) - ...
%                                                                                                              (sqrt((B*m)^2 - 1)/(B*m)))*((B*bodyDiameter/rootChord)^2)*acosh(1 + rootChord/(B*bodyDiameter)) - ...
%                                                                                                              ((B*m)/(1 + B*m))*acos(1/(B*m)));
% end
totalLift(1) = numFins*liftForce*KTB;
totalLift(2) = numFins*liftForce*KBT;

if attackAng==pi/180
    CLa(1) = CLaT*KTB;
    CLa(2) = CLaT*KBT;
else
    CLa = [0,0];
end

CPXt = pitchMoment/(liftForce*cos(attackAng));
centerPresX(1) = distNoseTip + CPXt;
centerPresX(2) = distNoseTip + rootChord/4 + ((sqrt(semispan^2 - bodyRadius^2)*acosh(semispan/bodyRadius) - semispan + pi*bodyRadius/2)/((bodyRadius/sqrt(semispan^2 - bodyRadius^2))*acosh(semispan/bodyRadius) + semispan/bodyRadius - pi/2) - bodyRadius)*(tan(leadSweep) - 0.25*(rootChord-tipChord)/semispan);
% MBT = (numFins*1.225*((343*Mach)^2)*attackAng*m*(rootChord^3)/(3*pi*B))*(sqrt(1 + 2*B*bodyDiameter/rootChord)*((2*m*B + 5)/(3*(m*B + 1)^2) + (B*bodyDiameter/rootChord)/(3*(B*m + 1)) - ((B*bodyDiameter/rootChord)^2)/(B*m)) + ...
%                                                                   (((m*B)^2 - 1)^-0.5)*((1 + B*bodyDiameter/rootChord)^3 - ((B*bodyDiameter/rootChord)^3)/((m*B)^2) - (1+m*B)^-2)*acos((1 + (B*bodyDiameter/rootChord)*(m*B + 1))/(m*B + (B*bodyDiameter/rootChord)*(m*B + 1))) + ...
%                                                                   (((B*bodyDiameter/rootChord)^3)/((m*B)^2))*acosh(1 + rootChord/(B*bodyDiameter)) - ((2*m*B + 5)/(3*(m*B + 1)^2)) - ((1 - (m*B + 1)^-2)/sqrt((m*B)^2 - 1))*acos(1/(m*B)));
%centerPresX(2) = distNoseTip + MBT/LBT;