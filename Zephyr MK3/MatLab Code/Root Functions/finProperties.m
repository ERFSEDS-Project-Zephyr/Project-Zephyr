function finProps = finProperties(finProps,finRes)
% Function that takes in the fin dimensions from the master file and
% computes additional properties useful for other calculation

% Sweep and Wedge Angles
finProps.leadSweepAngle = tan(finProps.sweepLength./finProps.height);
finProps.trailSweepAngle = tan((finProps.rootChord - finProps.tipChord - finProps.sweepLength)./finProps.height);
finProps.leadWedgeAngle = 2*tan(finProps.thickness./(2*finProps.leadTaperLength));
finProps.trailWedgeAngle = 2*tan(finProps.thickness./(2*finProps.trailTaperLength));

% Relevant Aerodynamic Areas
finProps.frontArea = finProps.thickness.*finProps.height;
finProps.sideArea = 0.5*finProps.height.*(finProps.rootChord + finProps.tipChord);
finProps.wettedArea = finProps.height.*(sqrt(4*finProps.leadTaperLength.^2 + finProps.thickness.^2)+...
                                        sqrt(4*finProps.trailTaperLength.^2 + finProps.thickness.^2)+...
                                        finProps.rootChord+finProps.tipChord-2*(finProps.leadTaperLength+finProps.trailTaperLength));

% Relevant Aerodynamic Quantities
finProps.aspectRatio = 8*(finProps.height.^2)./finProps.wettedArea;
finProps.meanAeroChord = (4/3)*(finProps.rootChord + finProps.tipChord - (finProps.rootChord.*finProps.tipChord)./(finProps.rootChord + finProps.tipChord));
finProps.midChordAngle = tan((finProps.rootChord - finProps.tipChord - 2*finProps.sweepLength)./(2*finProps.height));
finProps.midChordLength = sqrt((0.5*finProps.rootChord - 0.5*finProps.tipChord - finProps.sweepLength).^2 + finProps.height.^2);

% Center of Gravity Location
if ~exists('finRes','var')
    finRes = 100;
end
y = finProps.height.*linspace(0,1,finRes)';
dy = y(2)-y(1);

chordLength = finProps.tipChord + (finProps.rootChord-finProps.tipChord).*y./finProps.height;
V_total = 0; finProps.CG = zeros(2,size(finDims,2));
for k = 1:finRes
    V_section = dy*(chordLength(k,:) - 0.5*(finProps.leadTaperLength + finProps.trailTaperLength));
    V_total = V_total + V_section;

    finProps.CG(1,:) = finProps.CG(1,:) + V_section.*(chordLength(k,:).^2 - finProps.trailTaperLength.*chordLength(k,:) - (1/3)*finProps.leadTaperLength.^2)...
                                                   ./(2*chordLength(k,:) - finProps.leadTaperLength - finProps.trailTaperLength);
    finProps.CG(2,:) = finProps.CG(2,:) + V_section*(y(k) + 0.5*dy);
end
finProps.CG = finProps.CG./V_total;
finProps.CG(3,:) = (finProps.CG(1,:) - finProps.sweepLength.*finProps.CG(2,:)./finProps.height)...
                    ./(finProps.rootChord - (finProps.rootChord - finProps.tipChord).*finProps.CG(2,:)./finProps.height);