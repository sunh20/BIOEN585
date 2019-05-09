%% BIOEN585 Lab 6: Parameter Estimation
% Samantha Sun
% 20190508

% Measure strength of single receptor-ligand bonds under force

%% Question 1: Graph 4 sets of data
clear all; close all; clc;
load('dat.mat')

% plot
figure;
subplot(4,1,1)
bar(dat.bins,dat.f300)
ylabel('Occurences')
title('Bond Rupture Force @ 300 pN/s loading rate')

subplot(4,1,2)
bar(dat.bins,dat.f3000)
ylabel('Occurences')
title('Bond Rupture Force @ 3,000 pN/s loading rate')

subplot(4,1,3)
bar(dat.bins,dat.f30000)
ylabel('Occurences')
title('Bond Rupture Force @ 30,000 pN/s loading rate')

subplot(4,1,4)
bar(dat.bins,dat.ctrl)
xlabel('Rupture Force (pN)')
ylabel('Occurences')
title('Bond Rupture Force - Negative Control')

%% Question 2: Build ODE model

% guess parameters
k = 1;
r = 1;
fs = 1;
S_0 = 10;
guesses = [k,r,fs,S_0];

% fminsearch
[estimates, J] = fminsearch(@obj,dat,guesses);

% setup initial conditions + time
S_0 = 0; % don't know
tspan = [0 20];

% ode
[t,y] = ode45(@bondODE, tspan, S_0, [], params);











