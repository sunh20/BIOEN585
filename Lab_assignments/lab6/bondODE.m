% Describes the number of bonds surviving when the force = f

function dSdf = bondODE(f,S,params)
global dat

k = params(1);
fs = params(2);
a = params(3);

dSdf = zeros(4,1);
dSdf(1) = -k / dat.r(1) * exp(f/fs) * S(1);  % solves for S f300
dSdf(2) = -k / dat.r(2) * exp(f/fs) * S(2);  % solves for S f3000
dSdf(3) = -k / dat.r(3) * exp(f/fs) * S(3);  % solves for S f30000
dSdf(4) = -a * S(4) * exp(-a/f) / f^2;             % solves for N 

end