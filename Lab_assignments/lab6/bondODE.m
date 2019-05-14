% Describes the number of bonds surviving when the force = f
% Reminder that data is number of bonds that ruptured within range
% of forces

function dSdf = bondODE(f,y,params)
global dat
k = params(1);
fs = params(2);
r = params(3);

dSdf = -k / r * exp(f/fs) * y;

end