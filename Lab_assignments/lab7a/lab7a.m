%% Lab 7a: Modeling Brownian Dynamics
% Samantha Sun
% BIOEN 585
% 20190513

%% setup
clear all; close all; clc

% parameters
Dv = 2e-3;      % um^2/s
V_avg = 0.85;   % um/s
Lp = 111;       % um
params = [Dv,V_avg,Lp];

% time step
dt = 0.1;
tspan = 0:dt:1000;

%% Question 1: Simulate motion of single microtubule
% num tracks
N = 100;

tic
[X,Y,Q] = MT(params,tspan,dt,N);
toc

% plot
figure;
plot(X,Y)
xlabel('x position')
ylabel('y position')

% save('q1_dat.mat','X','Y','Q')

% part b: calculate V_avg and Dv from solved traces

% calculate contour length/distance traveled, v=d/t
S = sum(sqrt(diff(X).^2 + diff(Y).^2));
V_avg_est = mean(S./tspan(end));

% diffusion coefficient = dist diffused / (2t)
Dv_est = mean((S - mean(S)).^2./2/tspan(end));

% part c: calculate Lp - 
% avg of theta over all tracks at each time point vs
% expected value of s at each time point
%S_exp = cumsum(ones(length(tspan),1)*V_avg);
S_exp = V_avg*tspan;
Q_cosavg = mean(cos(Q),2);
Q_esLp = exp(-S_exp./Lp);

figure;
plot(S_exp,Q_cosavg,S_exp,Q_esLp)
xlabel('Expected contour length')
ylabel('some units')
legend('avg cos(theta)','Lp used')

Lp_est = nanmean(-S_exp./log(Q_cosavg)');

% part d: varying time steps
% time step
dt = [0.1, 0.5, 1, 2];
results = zeros(length(dt),3); % 3 params

for i = 1:length(dt)
    tspan = 0:dt(i):1000;

    % solve
    tic
    [X,Y,Q] = MT(params,tspan,dt(i),N);
    toc
    
    % save relevant values

    S = sum(sqrt(diff(X).^2 + diff(Y).^2));
    V_avg_est = mean(S./tspan(end));

    % diffusion coefficient = dist diffused / (2t)
    Dv_est = mean((S - mean(S)).^2./2/tspan(end));
    
    % Lp est
    S_exp = V_avg*tspan;
    Q_cosavg = mean(cos(Q),2);
    Lp_est = nanmean(real(-S_exp./log(Q_cosavg)'));
    
    results(i,:) = [V_avg_est, Dv_est, Lp_est];
end

%% Question 2: Simulate in microtubule
% run setup subsection

d = 5.5; % diameter of walls
N = 10;

tic
[X,Y,Q] = MT(params,tspan,dt,N,d);
toc

% plot
figure;
plot(X,Y)
hold on
plot(tspan,ones(size(tspan))*d/2,'k--','LineWidth',2)
plot(tspan,ones(size(tspan))*-d/2,'k--','LineWidth',2)
xlabel('x position')
ylabel('y position')
xlim([0 100])






