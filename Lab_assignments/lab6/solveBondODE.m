function rup = solveBondODE(guesses)
global dat
% determine the model behavior with the guesses values
fspan = 0:20:300;
S = zeros(length(fspan),length(dat.r));

% for each force condition, solve ODE
for i = 1:length(dat.r) 
    params = [guesses(1:2), dat.r(i)]
    S_0 = guesses(2+i);
    [f,S(:,i)] = ode23(@bondODE, fspan, S_0, [], params);
end

% bin number of ruptured bonds
rup = -diff(S); 
rup = rup(2:end,:);    % don't count 0-20 bin

% plot just f300
plot(fspan(2:end-1),rup(:,1),fspan(2:end-1),dat.f300,'o')
drawnow
xlabel('rupture force')
ylabel('occurences')

end