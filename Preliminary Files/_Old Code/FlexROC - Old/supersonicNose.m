function [CLa,centerPresX,liftForce] = supersonicNose(l,d,N,M)
rho = d/4 + (l^2)/d;
B = sqrt(M^2 - 1);
chi = [];

syms x r R(x);
R(x) = sqrt(rho^2 - (l-x)^2) + d/2 - rho;
dR = diff(R,x); e = dR(0);

t(1) = B*r/(x-chi(1));
C(1) = (e^2)/(sqrt(1 - (B*e)^2 + e^2)*asech(B*e));
phi(1) = -C(1)*x*(asech(t(1)) - sqrt(1 - t(1)^2));

for n = 2:N
    syms c;
    t(n) = B*r/(x-chi(n));
    phi(n) = -c*((x - chi(n))^2)*((1 + 0.5*t(n)^2)*asech(t(n)) - 1.5*sqrt(1 - t(n)^2));
    C(n) = solve(sum(diff(phi,r)) == dR(1 + sum(diff(phi,x))),c);
    phi(n) = -C(n)*((x - chi(n))^2)*((1 + 0.5*t(n)^2)*asech(t(n)) - 1.5*sqrt(1 - t(n)^2));
end