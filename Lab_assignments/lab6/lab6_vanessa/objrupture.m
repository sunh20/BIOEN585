function J = objrupture(guesses)

global data
global r

fspan = data(:,1);

guesses(1:2)=abs(guesses(1:2));

S0 = [guesses(3), guesses(4), guesses(5) guesses(6)];
options = [];

[f,s] = ode15s(@negODE, fspan, S0, options, guesses, r);

R = -diff(s);
R(:,1:3) = R(:,1:3) + R(:,4);
data_cut = data(1:end-1,:);

stdv_dat = sqrt(data_cut) + 1;

plot(f(1:end-1,1),R,data(:,1),data(:,2),'o',data(:,1),data(:,3),'*',data(:,1),data(:,4),'.',data(:,1),data(:,5),'s'); drawnow
% J = sum(((data_cut(:,2)-R(:,1))).^2)+sum(((data_cut(:,3)-R(:,2))).^2)+sum(((data_cut(:,4)-R(:,3))).^2)+sum(((data_cut(:,5)-R(:,4))).^2)
J = sum(((data_cut(:,2)-R(:,1))./stdv_dat(:,2)).^2)+sum(((data_cut(:,3)-R(:,2))./stdv_dat(:,3)).^2)+sum(((data_cut(:,4)-R(:,3))./stdv_dat(:,4)).^2)+sum(((data_cut(:,5)-R(:,4))./stdv_dat(:,5)).^2)
end