function [A,N,xcp,ycp] = noseAero(gamma,R,pinf,Tinf,M,alpha,d,l,t)

q = (gamma/2)*pinf*M^2;
B = sqrt(abs(1 - M^2));

f = l/d; F = 2*f + 1/(2*f);
area_ref = (pi/4)*d^2;
area = (pi/4)*d^2;
area_w = pi*l*(4*f - F + (F^2 / 4)*asin(2/F) + sqrt((F^2 / 4) - 1));
vol = area_w*t;

mu = (1.789*10^-5)*((Tinf/288.15)^(3/2))*((288.15 + 113))/(T + 113);
Re = (M*l/mu)*pinf*sqrt(gamma/(R*Tinf));
if M < 1
    A = ;
    N = (2*area/area_ref)*q*area/B;
elseif M > 1
    A = ;
    N = (2*area/area_ref)*q*area/B;
else
    A = ;
    N = (2*area/area_ref)*q*area/B;
end
xcp = l - vol/area;
ycp = 0;