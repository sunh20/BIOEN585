function J = objforce(guesses)

global data
global r

fspan = data(:,1);

guesses(1:2)=abs(guesses(1:2));

S0 = [guesses(3), guesses(4), guesses(5)];
options = [];
[f,s] = ode15s(@forceODE, fspan, S0, options, guesses, r);

R = -diff(s);
data_cut = data(1:end-1,:);
plot(f(1:end-1,1),R,data(:,1),data(:,2),'o',data(:,1),data(:,3),'*',data(:,1),data(:,4),'.'); drawnow
J = sum(((data_cut(:,2)-R(:,1))).^2)+sum(((data_cut(:,3)-R(:,2))).^2)+sum(((data_cut(:,4)-R(:,3))).^2);
end
