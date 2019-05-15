% Objective function (J)

function J = obj(guesses)
global dat
% solve ODE for all 4 conditions
%disp(guesses)
S_0 = guesses(4:end);
fspan = 20:20:300;
[f,S] = ode23s(@bondODE, fspan, S_0, [], guesses);
rup = -diff(S); 

% plot just f300
% plot(fspan(1:end-1),rup(:,1),fspan(1:end-1),dat.f300,'o',...
%     fspan(1:end-1),rup(:,2),fspan(1:end-1),dat.f3000,'o',...
%     fspan(1:end-1),rup(:,3),fspan(1:end-1),dat.f30000,'o')
% drawnow
% xlabel('rupture force')
% ylabel('occurences')
% legend('300','3000','30000')
% ylim([0 45])

% weighted least-squares
J = sum(((dat.f300-rup(:,1))).^2)/std(dat.f300) + ...
    sum(((dat.f3000-rup(:,2))).^2)/std(dat.f3000) + ...
    sum(((dat.f30000-rup(:,3))).^2)/std(dat.f30000); 
% + sum(((dat.ctrl-rup(:,4))).^2)/std(dat.ctrl);

if guesses(3) ~= 0
    J = J + sum(((dat.ctrl-rup(:,4))).^2)/std(dat.ctrl);
end