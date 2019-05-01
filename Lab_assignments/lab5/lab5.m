%% Lab 5: Stochastic chemical reaction system
% Samantha Sun
% BIOEN 585
% 20190501

%% set up system
clear all; close all; clc
% parameters
k1 = 1;
k2 = 1;
k3 = 100;
k4 = 1;
k5 = 100;
kD = 1;
params = [k1,k2,k3,k4,k5,kD];

% variables - IC all zero
A1_0 = 10;
A1_p_0 = 0;
A1_i_0 = 0;
A2_0 = 10;
E_0 = 0;
S_0 = 0;
IC = [A1_0, A1_p_0, A1_i_0,A2_0,E_0,S_0];

save("params.mat", 'params', 'IC')

%% Q1: Solve Artymov paper model
clear all; close all; clc

% import params and ICs
load("params.mat");

% % deterministic soln
% tspan = 0:0.01:20; 
% [t,y] = ode23s(@Tcellode,tspan,IC,[],params);

% Stochastic soln: exact
tspan = [0 20];
[t,y] = DSDEexact(@TCellRXN,tspan,IC,params);

% plot
figure;
subplot(2,1,1)
plot(t,y)
xlim([0 0.2])
ylim([0 11])
legend('A1','A1-prot','A1-inact','A2','E','S')
title('Fast-response')

subplot(2,1,2)
plot(t,y)
xlim([0 20])
ylim([0 11])
title('Slow-response')

%% Q2: Comparing solver efficiency, accuracy
clear all; close all; clc

% parameters
load("params.mat")
IC_set = [10, 30, 100, 300, 1000]; % for A1_0, A2_0

% efficiency: plot log time vs. A1_0 and A2_0 value
% 1. tau solver w/ relative tol = 0.1
% 2. tau solver w/ relative tol = 0.01
% 3. tau solver w/ relative tol = 0.001
% 4. exact solver

tspan = [0 20];
t_eff = zeros(length(IC_set),4); % num ICs x 4 solvers

for i = 1:length(IC_set)
    IC(1) = IC_set(i);
    IC(4) = IC_set(i);
    
    % tau solve @ 0.1
    tic
    [t,y] = DSDEtauleap(@TCellRXN,tspan,IC,0.1,params);
    t_eff(i,1) = toc;
    
    % tau solve @ 0.01
    tic
    [t,y] = DSDEtauleap(@TCellRXN,tspan,IC,0.01,params);
    t_eff(i,2) = toc;
    
    % tau solve @ 0.001
    tic
    [t,y] = DSDEtauleap(@TCellRXN,tspan,IC,0.001,params);
    t_eff(i,3) = toc;
    
    % exact
    tic
    [t,y] = DSDEleap(@TCellRXN,tspan,IC,params);
    t_eff(i,4) = toc;
end










