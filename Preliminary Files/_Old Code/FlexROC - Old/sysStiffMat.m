function Ksys = sysStiffMat(sectionDist,ID,OD,E,aeroIndex,liftForce,attackAng)
N = length(sectionDist);
I = pi*(OD^4 - ID^4)/64;
Kstr = zeros(2*N,2*N);
Kaero = zeros(2*N,2*N);

for i = 1:N-1
    l = sectionDist(i+1) - sectionDist(i);
    k = (E*I/l^3)*[12,6*l,-12,6*l;6*l,4*l^2,-6*l,2*l^2;-12,-6*l,12,-6*l;6*l,2*l^2,-6*l,4*l^2];
    L = [zeros(4,2*(i-1)),eye(4),zeros(4,2*(N-i-1))];
    Kstr = Kstr + L'*k*L;
end

if ~isempty(aeroIndex)
    for i = 1:length(aeroIndex)
        k = (liftForce(i)/attackAng)*[0,1;0,0];
        L = [zeros(2,2*(i-1)),eye(2),zeros(2,2*(N-i))];
        Kaero = Kaero + L'*k*L;
    end
end

Ksys = Kstr + Kaero;