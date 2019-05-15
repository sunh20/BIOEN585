%% Lab 7a: Modeling Brownian Dynamics
% Samantha Sun
% BIOEN 585
% 20190513

%% Question 1: Simulate motion of single microtubule
clear all; close all; clc

% parameters
Dv = 2e-3;      % um^2/s
V_avg = 0.85;   % um/s
Lp = 111;       % um
params = [Dv,V_avg,Lp];

% time step
dt = 0.1;
tspan = 0:dt:10;

% num tracks
N = 10;

[X,Y,Q] = MT(params,tspan,dt,N);



