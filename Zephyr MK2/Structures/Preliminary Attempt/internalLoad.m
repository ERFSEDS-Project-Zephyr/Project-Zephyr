function [N,V,M] = internalLoad(zcg,zcp,m,T,L,D,a)

g = 9.807;

W = m*g;

Fx = sum(L)*cos(a) + sum(D)*sin(a);
Fz = sum(L)*sin(a) - sum(D)*cos(a) - sum(W) + T;
M = sum(L.*zcp)*cos(a) + sum(D.*zcp)*sin(a);

Ns = ma - Fz; Vs = -Fx; Ms = M;

end

