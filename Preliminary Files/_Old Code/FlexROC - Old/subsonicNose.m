function [CLa,centerPresX] = subsonicNose(bodyOuterDiameter,noseLength)
bodyRadius = bodyOuterDiameter/2;
ogiveRadius = (noseLength^2 + bodyRadius^2)/(2*bodyRadius);
filledVolume = pi*((bodyRadius-ogiveRadius)*asin(noseLength/ogiveRadius)*ogiveRadius^2 + (bodyRadius^2 - 2*ogiveRadius*bodyRadius + 2*ogiveRadius^2)*noseLength + (bodyRadius-ogiveRadius)*noseLength*sqrt(ogiveRadius^2 - noseLength^2) - (noseLength^3)/3);
CLa = 2; % Eq. 3-66
centerPresX = noseLength - filledVolume/(pi*bodyRadius^2); % Eq. 3-89