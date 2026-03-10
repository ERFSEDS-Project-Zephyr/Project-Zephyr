function [A,N,xcp,ycp] = finAero(gamma,pinf,M,alpha,d,finDims,numSections)

q = (gamma/2)*pinf*M^2;
B = sqrt(abs(1 - M^2));

numFins = finDims(1); rootChord = finDims(2)/100; tipChord = finDims(3)/100; semispan = finDims(4)/100; finThickness = finDims(5)/100; leadTaperLength = finDims(6)/100; trailTaperLength = finDims(7)/100; leadEdgeLocation = finDims(8)/100;
leadSweepAng = finDims(9)*pi/180; leadTaperSweepAng = finDims(10)*pi/180; trailTaperSweepAng = finDims(11)*pi/180; trailSweepAng = finDims(12)*pi/180; leadWedgeAng = finDims(13)*pi/360; trailWedgeAng = finDims(14)*pi/360;

dy = semispan/numSections;
S = semispan*(rootChord+tipChord);
aspectRatio = (2*semispan)^2 / S;
meanChord = (2/3)*(rootChord + tipChord - (rootChord*tipChord)/(rootChord + tipChord));

CL = 0; CN = 0;
CD = 0; CA = 0;
if M<1
    for y = 0:dy:semispan
        chord = rootChord + (tipChord - rootChord)*y/semispan;
        L = L + (2*pi*alpha / B)*chord*q;
        D = (0.3*finThickness/chord + (2*pi*alpha/B)^2 / (pi*aspectRatio))*chord*q;
    end

    A = D*cos(alpha) - L*sin(alpha);
    N = D*sin(alpha) + L*cos(alpha);

    ycp = d/2 + (semispan/3)*(1 + tipChord/(rootChord + tipChord));
    xcp = leadEdgeLocation + (ycp - d/2)*tan(leadSweepAng) + meanChord/4;
elseif M>1
    
else
    A = Inf;
    N = Inf;
    xcp = Inf;
    ycp = Inf;
end