% Objective function (J)

function J = obj(dat,guesses)

% determine the model behavior with the guesses values:
tspan = 0:0.1:10;
S_0 = guesses(4);
y = cell(4,1);

for i = 1:length(dat.f)
    params = [guesses(1:3), dat.f(i)];
    [t,y{i}] = ode45(@bondODE, tspan, S_0, [], params);
end

% weighted least-squares hi this is wrong i stopped here
J = sum(((dat.f300-y{1})).^2) + sum(((dat.f3000-y{2})).^2) + ...
    sum(((dat.f30000-y{3})).^2) + sum(((dat.ctrl-y{4})).^2);
