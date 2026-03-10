function [CLa,centerPresX] = subsonicFin(bodyRadius,finDims,Mach)

numFins = finDims(1); rootChord = finDims(2)/100; 
tipChord = finDims(3)/100; semispan = finDims(4)/100; 
distNoseTip = finDims(7)/100; leadSweep = finDims(8)*pi/180;  

t = (semispan + bodyRadius)/(bodyRadius);
KTB = (2/(pi*(1 - t^-1)^2))*((1 + t^-4)*(0.5*atan(0.5*(t - t^-1))) - (t^-2)*(t - t^-1 + 2*atan(1/t)));
KBT = ((1 + 1/t)^2) - KTB;

midSweep = atan(tan(leadSweep) - 0.5*(rootChord - tipChord)/semispan);

finArea = 0.5*(rootChord + tipChord) * semispan;
AspectRatio = (semispan)^2 / finArea;

refArea = pi*(bodyRadius)^2;
CLaT = (numFins* pi * (AspectRatio * finArea / refArea)) / (2 + sqrt(4 + ( (sqrt(abs(1 - Mach^2))) * AspectRatio / cos(midSweep) )^2 ));

xT = semispan * tan(leadSweep);
xBarT = distNoseTip + (xT/3) * ((rootChord + 2 * tipChord) / (rootChord + tipChord)) + (1/6) * (rootChord + tipChord - (rootChord * tipChord) / (rootChord + tipChord)) ;

CLa(1) = CLaT * KTB;
CLa(2) = CLaT * KBT;

centerPresX(1) = xBarT;
centerPresX(2) = distNoseTip + rootChord/4 + ((sqrt(semispan^2 - bodyRadius^2)*acosh(semispan/bodyRadius) - semispan + pi*bodyRadius/2)/((bodyRadius/sqrt(semispan^2 - bodyRadius^2))*acosh(semispan/bodyRadius) + semispan/bodyRadius - pi/2) - bodyRadius)*(tan(leadSweep) - 0.25*(rootChord-tipChord)/semispan);