function dsdf = negODE(f,S,est,r)

k = est(1);
fs = est(2);
a = est(7);

dsdf = zeros(4,1);

dsdf(1) = -(k/r(1)) * exp(f/fs) * S(1);
dsdf(2) = -(k/r(2)) * exp(f/fs) * S(2);
dsdf(3) = -(k/r(3)) * exp(f/fs) * S(3);
dsdf(4) = -a * S(4);
end
