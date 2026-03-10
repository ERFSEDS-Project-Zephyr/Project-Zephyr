function [CLa,centerPresX,liftForce,stability] = subsonicRocket(bodyOuterDiameter,noseLength,numFinSets,Mach,attackAng,finDims,rocketCG)

[CLa(1),centerPresX(1)] = subsonicNose(bodyOuterDiameter,noseLength);

for i=1:numFinSets
    [CLa(2*i:2*i + 1),centerPresX(2*i:2*i + 1)] = subsonicFin(bodyOuterDiameter/2,finDims(:,i),Mach);
end
liftForce = numFinSets*(0.5*1.225*(Mach*340.3))*CLa*attackAng;

rocketCenterPres = sum(centerPresX.*CLa)/sum(CLa);
stability = (rocketCenterPres - rocketCG)/bodyOuterDiameter;