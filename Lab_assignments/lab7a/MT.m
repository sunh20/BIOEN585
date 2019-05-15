function [X,Y,Q] = MT(params,tspan,dt,N)
% solves for N number of microtubule tracks over time t

% define needed terms
Dv = params(1);
V_avg = params(2);
Lp = params(3);

R_std = sqrt(2*Dv*dt);
Q_std = sqrt(V_avg*dt/Lp);

% variables - IC = 0
X = zeros(length(tspan),length(N));
Y = zeros(length(tspan),length(N));
Q = zeros(length(tspan),length(N));

% calculate tracks
for t = 2:length(tspan)
    dQ = normrnd(0,Q_std,[N,1]);
    r = normrnd(V_avg*dt,R_std,[N,1]);
    
    Q(t,:) = Q(t-1,:) + dQ;     % new sigma term
    dx = r*cos(Q(t,:));
    dy = r*sin(Q(t,:));
    X(t,:) = X(t-1,:) + dx;
    Y(t,:) = Y(t-1,:) + dy;
    
end


end