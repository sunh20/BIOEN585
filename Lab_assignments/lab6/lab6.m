%% BIOEN585 Lab 6: Parameter Estimation
% Samantha Sun
% 20190508

% Measure strength of single receptor-ligand bonds under force

%% setup if needed
if ~isfile('dat.mat')
    disp('Setting up data...')
    setup
end
%% Question 1: Graph 4 sets of data
clear all; close all; clc;
load('dat.mat')
global dat

plotCompare([],0)

%% Question 2: Build ODE model
% Question 3: solve ODE model with initial guesses

% guess parameters
k = 30; % s^-1
fs = 30; % pN
S_0 = [sum(dat.f300),sum(dat.f3000),sum(dat.f30000),sum(dat.ctrl)];

guesses = [k,fs,S_0];

% solve ODE using guess
rup = solveBondODE(guesses);

% plot
plotCompare(rup)

%% Question 4: Estimate parameters (w/o neg control)

% use fminsearch to find parameters
tic
[estimates, J] = fminsearch(@obj,guesses);
toc
disp(estimates)

% solve ODE using estimated parameters


%[t,y] = ode45(@bondODE, tspan, S_0, [], params);











