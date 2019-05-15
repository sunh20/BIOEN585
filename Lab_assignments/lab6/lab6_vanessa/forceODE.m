function dsdf = forceODE(f,S,est,r)

k = est(1);
fs = est(2);

dsdf = zeros(3,1);

dsdf(1) = -(k/r(1)) * exp(f/fs) * S(1);
dsdf(2) = -(k/r(2)) * exp(f/fs) * S(2);
dsdf(3) = -(k/r(3)) * exp(f/fs) * S(3);
end