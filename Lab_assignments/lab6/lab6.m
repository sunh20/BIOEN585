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
%plotCompare([],0)

%% Question 2: Build ODE model
% Question 3: solve ODE model with initial guesses

% guess parameters + IC
k = 20; % s^-1
fs = 30; % pN
S_0 = [sum(dat.f300);sum(dat.f3000);sum(dat.f30000);sum(dat.ctrl)];
params = [k,fs];

% solve ODE using guess
fspan = 20:20:300;
[f,S] = ode23s(@bondODE, fspan, S_0, [], params);

% num ruptured
rup = -diff(S); 

% plot
plotCompare(rup)

%% Question 4: Estimate parameters (w/o neg control)

% guess params
k = 0.1; % s^-1
fs = 10; % pN
S_0 = [80,80,80,80];

guesses = [k,fs,0,S_0]; % 0 padding for later variable

% use fminsearch to find parameters
tic
[estimates, J] = fminsearch(@obj,guesses);
toc
disp(estimates)

% solve ODE using estimated parameters
fspan = 20:20:300;
params = [estimates(1),estimates(2),0];
S_0 = estimates(4:end);
[f,S] = ode23s(@bondODE, fspan, S_0, [], params);

% num ruptured
rup = -diff(S); 

% plot
plotCompare(rup)

% save
save('estimates.mat','estimates','J')

%% Question 5: adding negative control

clear all; close all; clc
load('dat.mat')
global dat

% guess parameters
k = 0.1; % s^-1
fs = 10; % pN
S_0 = [80,80,80,80];
a = 0.1;

guesses = [k,fs,a,S_0];

% use fminsearch to find parameters
tic
[estimates, J] = fminsearch(@obj,guesses);
toc
disp(estimates)

% solve ODE using estimated parameters
fspan = 20:20:300;
params = [estimates(1),estimates(2),estimates(3)];
S_0 = estimates(4:end);
[f,S] = ode23s(@bondODE, fspan, S_0, [], params);

% num ruptured
rup(:,4) = -diff(S(:,4));
rup(:,1:3) = -diff(S(:,1:3) + S(:,4)); 

% plot
plotCompare(rup)

% save
save('estimates_neg.mat','estimates','J')


