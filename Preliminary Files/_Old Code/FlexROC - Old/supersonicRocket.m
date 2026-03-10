function [CLa,centerPresX,liftForce,stability] = supersonicRocket(bodyOuterDiameter,noseLength,numFinSets,finResolution,Mach,attackAng,finDims,rocketCG)

[CLa(1),centerPresX(1)] = subsonicNose(bodyOuterDiameter,noseLength);
for i=1:numFinSets
    [CLa(2*i:2*i + 1),~,~] = supersonicFin(bodyOuterDiameter,finResolution,Mach,pi/180,finDims(:,i),0);
    [~,centerPresX(2*i:2*i + 1),liftForce(2*i:2*i+1)] = supersonicFin(bodyOuterDiameter,finResolution,Mach,attackAng,finDims(:,i),0);
end
rocketCenterPres = sum(centerPresX.*CLa)/sum(CLa);
stability = (rocketCenterPres-rocketCG)/bodyOuterDiameter;