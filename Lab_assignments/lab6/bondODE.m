% Describes the number of bonds surviving when the force = f
% Reminder that data is number of bonds that ruptured within range
% of forces

function dSdf = bondODE(f,y,params)
global dat
k = params(1);
fs = params(2);
r = params(3);
neg = params(4); % tells us if negative case
a = params(5);

if ~neg
    dSdf = -k / r * exp(f/fs) * y;  % solves for S (y)
else
    dSdf = -a * y * exp(-a/f);      % solves for N (y)
end

end