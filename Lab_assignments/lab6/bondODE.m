% Describes the number of bonds surviving when the force = f
% Reminder that data is number of bonds that ruptured within range
% of forces

function dSdf = bondODE(t,y,params)

k = params(1);
r = params(2);
f = params(3);
fs = params(4);

dSdf = -k / r * exp(f/fs) * y;

end