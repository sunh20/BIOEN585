% Objective function (J)

function J = obj(guesses)
global dat
% solve ODE for all 4 conditions
tic
disp(guesses)
rup = solveBondODE(guesses,0);
toc

% weighted least-squares
J = sum(((dat.f300-rup(:,1))).^2)/std(dat.f300) + ...
    sum(((dat.f3000-rup(:,2))).^2)/std(dat.f3000) + ...
    sum(((dat.f30000-rup(:,3))).^2)/std(dat.f30000); 
% + sum(((dat.ctrl-rup(:,4))).^2)/std(dat.ctrl);
