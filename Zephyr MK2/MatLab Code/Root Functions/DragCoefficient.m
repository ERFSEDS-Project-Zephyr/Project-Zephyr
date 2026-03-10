function CD = DragCoefficient_GH(bodyDims,finDims,massDot,Mach,gamma,R,T,rho)

addpath('..')
bodyOuterDiameter = bodyDims(1)/100;
bodyInnerDiameter = bodyDims(2)/100;
noseLength = bodyDims(3)/100;
bodyDensity = bodyDims(4)*1000;
finResolution = bodyDims(5);
numStages = bodyDims(6);
bodyLength = sum(bodyDims(7:6+numStages))/100;
rocketLength = noseLength + bodyLength;

numFins = finDims(1,:);
rootChord = finDims(2,:)/100;
tipChord = finDims(3,:)/100;
semispan = finDims(4,:)/100;
thickness = finDims(5,:)/100;
leadSweep = finDims(9,:)*pi/180;
trailSweep = finDims(12,:)*pi/180;
leadWedge = finDims(13,:)*pi/360;
del_LE = pi/2 - leadWedge;

MAC = rootChord + 0.5*(tan(trailSweep) - tan(leadSweep)).*semispan + (((tan(trailSweep) - tan(leadSweep)).*semispan).^2)./(12*(rootChord + 0.5*(tan(trailSweep) - tan(leadSweep)).*semispan));
Sf = semispan.*(rootChord+tipChord);
Ss = Sf + thickness.*tipChord + thickness.*semispan./cos(leadSweep) + thickness.*semispan./cos(trailSweep);

V = Mach*sqrt(gamma*R*T);
q = 0.5*rho*V^2;

if massDot == 0 && Mach < 1
    CD_BB = 0.12 + 0.13*Mach^2;
    CD_BF = 0.053*(rocketLength/bodyOuterDiameter)*(Mach/(q*rocketLength))^0.2;
    CD_BT = CD_BB + CD_BF;

    CD_FF = numFins .* (0.0133 * (Mach./(q * MAC)).^0.2).* (2 * Sf ./ Ss);
    CD_FT = sum(CD_FF);

    CD = CD_BT + CD_FT;
elseif massDot == 0 && Mach > 1
    CD_BW = (1.586 + 1.834/Mach^2)*(atan(0.5/(noseLength/bodyOuterDiameter)))^1.69;
    CD_BB = 0.25/Mach;         
    CD_BF = 0.053*(rocketLength/bodyOuterDiameter)*(Mach/(q*rocketLength))^0.2;
    CD_BT = CD_BW + CD_BB + CD_BF;

    M_del = Mach .* cos(leadSweep);
    CD_FF = numFins .* (0.0133 * (Mach./(q * MAC)).^0.2) .* (2 * Sf ./ Ss);
    CD_FW = numFins .* (1.429./M_del.^2) .* (((1.2 * M_del.^2).^3.5) .* (2.4 ./ (2.8 * M_del.^2 - 0.4)).^2.5 - 1)...
        .* sin(del_LE).^2 .* cos(leadSweep) .* thickness .* semispan ./ Ss;
    CD_FT = sum(CD_FF) + sum(CD_FW);

    CD = CD_BT + CD_FT;
elseif massDot ~= 0 && Mach < 1
    CD_BF = 0.053*(rocketLength/bodyOuterDiameter)*(Mach/(q*rocketLength))^0.2;
    CD_FF = sum(numFins .* (0.0133 * (Mach./(q * MAC)).^0.2) .* (2 * Sf ./ Ss));
    CD = CD_BF + CD_FF;
else
    CD_BW = (1.586 + 1.834/Mach^2)*(atan(0.5/(noseLength/bodyOuterDiameter)))^1.69;
    CD_BF = 0.053*(rocketLength/bodyOuterDiameter)*(Mach/(q*rocketLength))^0.2;
    CD_FF = numFins .* (0.0133 * (Mach./(q * MAC)).^0.2) .* (2 * Sf ./ Ss);
    M_del = Mach .* cos(leadSweep);
    CD_FW = numFins .* (1.429./M_del.^2) .* (((1.2 * M_del.^2).^3.5) .* (2.4 ./ (2.8 * M_del.^2 - 0.4)).^2.5 - 1)...
        .* sin(del_LE).^2 .* cos(leadSweep) .* thickness .* semispan ./ Ss;
    CD = CD_BW + CD_BF + sum(CD_FF) + sum(CD_FW);
end