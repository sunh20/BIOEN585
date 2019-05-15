function rup = solveBondODE(guesses,neg)

global dat
% determine the model behavior with the guesses values
fspan = 0:20:300;
S = zeros(length(fspan),length(dat.r));
S_0_ref = [sum(dat.f300),sum(dat.f3000),sum(dat.f30000),sum(dat.ctrl)];

% for each force condition, solve ODE
for i = 1:length(dat.r)-1
    params = [guesses(1:2), dat.r(i),0,0,0];
    S_0 = S_0_ref(i);
    [f,S(:,i)] = ode23(@bondODE, fspan, S_0, [], params);
end

if ~neg % not solving for negative case
    % bin number of ruptured bonds
    rup = -diff(S); 
    rup = rup(2:end,:);    % don't count 0-20 bin
else    % with negative case
    for i = length(dat.r)
        params = [guesses(1:2), dat.r(i),1,guesses(4)];
        N_0 = guesses(3);
        [f,S(:,i)] = ode23(@bondODE, fspan, N_0, [], params);
    end

    % bin number of ruptured bonds
    rup(:,4) = -diff(S(:,4));
    rup(:,1:3) = -diff(S(:,1:3) + S(:,4)); 
    rup = rup(2:end,:);    % don't count 0-20 bin
end


% plot just f300
plot(fspan(2:end-1),rup(:,1),fspan(2:end-1),dat.f300,'o')
drawnow
xlabel('rupture force')
ylabel('occurences')

end