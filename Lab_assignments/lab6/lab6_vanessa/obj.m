function J = obj(guesses)
global data
% determine the model behavior with the guesses values:
tspan = data(:,1)-1900;
y0 = [guesses(5) guesses(6)];
options=[];
[t,y] = ode23s(@lotkavolterraODE,tspan,y0,options,guesses);
plot(t,data(:,2),'o',t,y(:,1),t,data(:,3),'*',t,y(:,2)) ; drawnow
% determine the weighted least-squares by comparing model to data:
J = sum(((data(:,2)-y(:,1))).^2)+sum(((data(:,3)-y(:,2))).^2);
% note that the above function assumes a certain error model. 
