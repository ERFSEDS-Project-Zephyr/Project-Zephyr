function M2 = prandtlMeyer(M1,del,g)
syms M
eqn = del + sqrt((g+1)/(g-1))*atan(sqrt(M1^2 - 1)/sqrt((g+1)/(g-1))) - atan(sqrt(M1^2 - 1))...
         == sqrt((g+1)/(g-1))*atan(sqrt(M^2 - 1)/sqrt((g+1)/(g-1))) - atan(sqrt(M^2 - 1));
M2 = abs(double(solve(eqn,M)));