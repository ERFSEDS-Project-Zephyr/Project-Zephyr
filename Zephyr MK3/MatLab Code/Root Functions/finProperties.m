function finProps = finProperties(finProps)
% Function that takes in the fin dimensions from the master file and
% computes additional properties useful for other calculation

% Sweep and Wedge Angles
finProps.leadSweepAngle = atand(finProps.sweepLength./finProps.height);
finProps.trailSweepAngle = atand((finProps.rootChord - finProps.tipChord - finProps.sweepLength)./finProps.height);
finProps.leadWedgeAngle = 2*atand(finProps.thickness./(2*finProps.leadTaperLength));
finProps.trailWedgeAngle = 2*atand(finProps.thickness./(2*finProps.trailTaperLength));

% Relevant Aerodynamic Areas
finProps.frontArea = finProps.thickness.*finProps.height;
finProps.sideArea = 0.5*finProps.height.*(finProps.rootChord + finProps.tipChord);
finProps.wettedArea = finProps.height.*(sqrt(4*finProps.leadTaperLength.^2 + finProps.thickness.^2)+...
                                        sqrt(4*finProps.trailTaperLength.^2 + finProps.thickness.^2)+...
                                        finProps.rootChord+finProps.tipChord-2*(finProps.leadTaperLength+finProps.trailTaperLength));

% Relevant Aerodynamic Quantities
finProps.aspectRatio = 8*(finProps.height.^2)./finProps.sideArea;
finProps.meanAeroChord = (4/3)*(finProps.rootChord + finProps.tipChord - (finProps.rootChord.*finProps.tipChord)./(finProps.rootChord + finProps.tipChord));
finProps.midChordAngle = atand((2*finProps.sweepLength + finProps.tipChord - finProps.rootChord)./(2*finProps.height));
finProps.midChordLength = sqrt((0.5*finProps.rootChord - 0.5*finProps.tipChord - finProps.sweepLength).^2 + finProps.height.^2);

% Center of Gravity Location
taperRatio = finProps.tipChord./finProps.rootChord;
sweepRatio = finProps.sweepLength./finProps.rootChord;
leadTaper = finProps.leadTaperLength./finProps.rootChord;
trailTaper = finProps.trailTaperLength./finProps.rootChord;

finProps.CG(2,:) = (2*(1 + 2*taperRatio) - 3*(leadTaper + trailTaper))./(6*(1 + taperRatio - leadTaper - trailTaper));
finProps.CG(3,:) = 0.5*(1 - (leadTaper.^2 - trailTaper.^2 - 1.5*(1 + taperRatio).*leadTaper)./(1 + taperRatio + taperRatio.^2 - 1.5*(1 + taperRatio).*(leadTaper + trailTaper)));
finProps.CG(1,:) = finProps.CG(3,:) + (sweepRatio - (1 - taperRatio).*finProps.CG(3,:)).*finProps.CG(2,:);