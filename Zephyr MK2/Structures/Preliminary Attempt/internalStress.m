function [normStress,shearStress,normFail,shearFail] = internalStress(nz,nx,alpha,l,do,di,a,T,zcg,m,zcp,L,D)
zcp = l - zcp; I = (pi/64)*(do^4 - di^4); E = 11*10^9;
maxNorm = max(0.25*(pi^2)*E*I/(l^2),380*10^6); maxShear = 380*10^6;

W = 9.807*m;
Ns = sum(L*sin(alpha)) - sum(D*cos(alpha)) + sum(W) + T - sum(m*a); % sum(Fz) = Lz - Dz - W + T - Ns = ma
Vs = sum(L*cos(alpha)) + sum(D*sin(alpha)); % sum(Fx) = Lx + Dx - Vs = 0
Ms = sum(L*cos(alpha).*zcp) + sum(D*sin(alpha).*zcp); % sum(M) = -Lx*zcp - Dx*zcp + Ms = 0

j = 1;
for z = linspace(l/nz,l,nz)
    i_cg = zcg<=z; i_cp = zcp<=z;
    
    N(j) = sum(m*a) - sum(L(i_cp)*sin(alpha)) + sum(D(i_cp)*cos(alpha)) + sum(W(i_cg)) - T + Ns; % Lz - Dz - W + T - Ns + N = ma
    V(j) = Vs - sum(L(i_cp)*cos(alpha)) - sum(D(i_cp)*sin(alpha)); % Lx + Dx - Vs + V = 0
    M(j) = Ms - V(j)*z - sum(L(i_cp)*cos(alpha).*zcp(i_cp)) - sum(D(i_cp)*sin(alpha).*zcp(i_cp)); % -Lx*zcp - Dx*zcp + Ms - V*x - M = 0
    
    k = 1;
    for x = -linspace(-do/2,do/2,nx)
        r = abs(x); Ro = do/2; Ri = di/2;
        
        Q = ((2/3)*(Ro^2 - r^2)^(3/2)) - ((2/3)*(Ri^2 - r^2)^(3/2));
        t = 2*Ro*sqrt(1 - (r/Ro)^2) - 2*Ri*sqrt(1 - (r/Ri)^2);
        
        normStress(j,k) = N(j)/((pi/4)*(do^2 - di^2)) - M(j)*x/I;
        shearStress(j,k) = V(j)*Q/(I*t);
        
        k = k + 1;
    end
    j = j + 1;
end

normFail = abs(normStress)>=maxNorm;
shearFail = abs(shearStress)>=maxShear;