function [CLa,centerPresX,liftForce,dragForce,stability] = subsonicRocket(bodyOuterDiameter,noseLength,numFinSets,Mach,attackAng,finDims,rocketCG)

[CLa(1),centerPresX(1),liftForce(1),dragForce(1)] = subsonicNose(bodyOuterDiameter,noseLength,Mach,attackAng);

for i=1:numFinSets
    [CLa(2*i:2*i + 1),centerPresX(2*i:2*i + 1),liftForce(2*i:2*i+1),dragForce(i+1)] = subsonicFin(bodyOuterDiameter/2,finDims(:,i),Mach,attackAng);
end

rocketCenterPres = sum(centerPresX.*CLa)/sum(CLa);
stability = (rocketCenterPres - rocketCG)/bodyOuterDiameter;