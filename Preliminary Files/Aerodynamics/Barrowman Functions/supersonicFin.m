function [CLa,centerPresX,totalLift,totalDrag] = supersonicFin(bodyDiameter,finResolution,Mach,attackAng,finDims,k)
liftForce = 0;
dragForce = 0;
pitchMoment = 0;

numFins = finDims(1); rootChord = finDims(2)/100; tipChord = finDims(3)/100; semispan = finDims(4)/100; finThickness = finDims(5)/100; leadTaper = finDims(6)/100; trailTaper = finDims(7)/100; distNoseTip = finDims(8)/100;
leadSweep = finDims(9)*pi/180; leadTaperSweep = finDims(10)*pi/180; trailTaperSweep = finDims(11)*pi/180; trailSweep = finDims(12)*pi/180; leadWedge = finDims(13)*pi/360; trailWedge = finDims(14)*pi/360;

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
    regionNormalLength(1) = sectionLeadLength*sin(incl(1));
    regionNormalLength(2) = sectionMidLength*sin(incl(2));
    regionNormalLength(3) = sectionTrailLength*sin(incl(3));
    regionNormalLength(4) = sectionLeadLength*sin(incl(4));
    regionNormalLength(5) = sectionMidLength*sin(incl(5));
    regionNormalLength(6) = sectionTrailLength*sin(incl(6));

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
        regionDragForce(I) = regionPresCoeff(I)*regionNormalLength(I)*(1 - 0.5*machConeRatio(I));
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
    sectionDragForce = 0;
    for I=1:6
        sectionLiftForce = sectionLiftForce + regionLiftForce(I);
        sectionDragForce = sectionDragForce + regionDragForce(I);
    end
    liftForce = liftForce + sectionLiftForce;
    dragForce = dragForce + sectionDragForce;
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

% KBT = ((8*B*bodyDiameter/(pi*rootChord*sqrt((m*B)^2 - 1))) / (B*CLaT*(tipChord/rootChord + 1)*(s/bodyRadius - 1)))*...
%       (((1 + m*rootChord/bodyDiameter)^2)*acos((m*B + rootChord/(B*bodyDiameter))/(1 + m*rootChord/bodyDiameter)) -...
%       ((m*B*rootChord/(B*bodyDiameter))^2)*acos((m*B)^-1) +...
%       (m*B*(rootChord/(B*bodyDiameter))^2)*sqrt((m*B)^2 - 1)*asin(B*bodyDiameter/rootChord) -...
%       sqrt((m*B)^2 - 1)*acosh(rootChord/(rootChord/(B*bodyDiameter))));

totalLift(1) = liftForce*KTB*refArea;
totalLift(2) = liftForce*KBT*refArea;
totalDrag(1) = dragForce*refArea;
totalDrag(2) = dragForce*refArea;

if attackAng==pi/180
    CLa(1) = CLaT*KTB;
    CLa(2) = CLaT*KBT;
else
    CLa = [0,0];
end

CPXt = pitchMoment/(liftForce*cos(attackAng));
centerPresX(1) = distNoseTip + CPXt;
centerPresX(2) = distNoseTip + rootChord/4 + ((sqrt(semispan^2 - bodyRadius^2)*acosh(semispan/bodyRadius) - semispan + pi*bodyRadius/2)/((bodyRadius/sqrt(semispan^2 - bodyRadius^2))*acosh(semispan/bodyRadius) + semispan/bodyRadius - pi/2) - bodyRadius)*(tan(leadSweep) - 0.25*(rootChord-tipChord)/semispan);
% q = (1/2)*(1.225)*(343*Mach)^2;
% 
% bdcr = B * bodyDiameter / rootChord; mb = m * B; mb1 = mb + 1; mb_inv = mb^-1; mbsqrt = sqrt(mb^2 - 1) ^ -1; mb23 = (2 * mb1 + 3) / (3 * mb1^2);
% 
% front = (4 * q * attackAng * m * rootChord^3) / (3 * pi * B);
% term1 = sqrt(1 + 2 * bdcr) * (mb23 + bdcr/(3 * mb1) - (bdcr^2)*mb_inv);
% term2 = mbsqrt * ...
%         ((1 + bdcr)^3 - (bdcr^3 * mb_inv^2) - (mb1^-2)) * ...
%         acos((1 + bdcr * mb1)/(mb + bdcr * mb1));
% term3 = bdcr^3 * mb_inv^2 * acosh(1 + bdcr^-1);
% term4 = mb23 + (1 - mb1^-2) * mbsqrt * acos(mb_inv);
% 
% MBT = front * ( term1 + term2 + term3 - term4 );
% centerPresX(2) = MBT/(KBT*liftForce);