function Msys = sysMassMat(sectionDist,sectionPointMass,ID,OD,rho)
N = length(sectionDist);
A = pi*(OD^2 - ID^2)/4;
Mlen = zeros(2*N,2*N);
Mpoint = zeros(2*N,2*N);

for i = 1:N-1
    l = sectionDist(i+1) - sectionDist(i);
    m = (rho*A*l/420)*[156 22*l 54 -13*l; 22*l 4*l^2 13*l -3*l^2; 54 13*l 156 -22*l; -13*l -3*l^2 -22*l 4*l^2];
    L = [zeros(4,2*(i-1)),eye(4),zeros(4,2*(N-i-1))];
    Mlen = Mlen + L'*m*L;
end

for i = 1:N
    m = sectionPointMass(i)/2;
    if i == 1
        L = [1,zeros(1,2*N-1)];
    elseif i == N
        L = [zeros(1,2*N-1),1];
    else
        L = [zeros(2,2*(i-1)),eye(2),zeros(2,2*(N-i))];
    end
    Mpoint = Mpoint + L'*m*L;
end

Msys = Mlen + Mpoint;