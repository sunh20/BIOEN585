%% Lab 3: Non-linear systems
% Samantha Sun
% BIOEN 585
% 20190417

clear all; close all; clc
plotOn = 0;

%% Question 1: Stability Analysis

% solve for nullclines + equilibrium points - plot on phase plane plot
a1 = 3.5;
a2 = 4;
beta = 2;

% polynomial eqn - solving for u, v
u_co = [2, -7, 4, -14, 34, -7];   % polynomial coefficients (found analytically)
null_u = roots(u_co);
null_u = null_u (null_u > 0);
null_v = a2 ./ (1 + null_u.^2);

% plot phase plane
u_range = 0:0.01:10;
v_range = 0:0.01:10;

plot(null_u,null_v,'--o')
xlabel('u')
ylabel('v')
xlim([0 5])
ylim([0 5])
title('Phase Plane Plot of u, v, equilibrium points')