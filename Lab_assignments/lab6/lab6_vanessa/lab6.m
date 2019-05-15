%% BIOEN 485 Lab 6: Parameter Estimation

%% 1. Graph the four sets of data
clear all; close all; clc;
global data
data = ...
[
20 5 8 4 9
40 4 7 5 6
60 5 8 4 6
80 5 5 4 4
100 5 3 6 3
120 14 7 5 3
140 25 11 8 3
160 9 27 8 3
180 2 9 20 2
200 2 3 11 1
220 1 2 4 1
240 1 0 1 2
260 1 0 2 1
280 0 0 2 1
];

mid_bin = data(:,1) + 10;                                                   %create middle of force bin vector
titles = ["300 pN/s", "3000 pN/s", "30,000 pN/s", "negative control"];      %Vector of strings for title

for i=1:4
    subplot(2,2,i)
    bar(mid_bin,data(:,i+1));
    title(titles(i))
    set(gca,'fontsize',14)
end

%% 3. Solve the model with initial guesses

global r
r = [300 3000 30000];                                                       %Vector of r values

IC = zeros(1,4);                                                            %Initial conditions equal sum of data in that category
for i=1:4
    IC(i) = sum(data(:,i+1));
end

guesses = [20 30 IC(1) IC(2) IC(3) IC(4) 0.01];                             %Set up vector for parameter and initial condition guesses

options = [];

fspan = [0 300];
[f,s] = ode23s(@forceODE, fspan, IC(1:3), options, guesses, r);             %Evaluate ODE with initial guesses

R = -diff(s);                                                               %Calculate R and plot

plot(f(2:end,1),R,mid_bin,data(:,2),'o',mid_bin,data(:,3),'*',mid_bin,data(:,4),'.')
xlabel('Bins')
ylabel('Ruptures')
legend('300 pN/s sim','3000 pN/s sim','30000 pN/s sim','300 pN/s data','3000 pN/s data','30000 pN/s data')

%% 4a. Estimate the parameters ignoring the negative control

guesses(1) = 0.1;                                                           %Use closer initial cnditions
guesses(2) = 11;

[estimates, J] = fminsearch(@objforce, guesses);                            %Minimize objective function
disp(estimates);

fspan = data(:,1);

S0 = [estimates(3), estimates(4), estimates(5)];

[f,s] = ode23s(@forceODE, fspan, S0, options, estimates, r);                %Evaluate ODE with minimized guesses

R = -diff(s);                                                               %Evaluate R and plot

plot(f(1:end-1,1),R,data(:,1),data(:,2),'o',data(:,1),data(:,3),'*',data(:,1),data(:,4),'.');
legend('300 pN/s sim','3000 pN/s sim','30000 pN/s sim','300 pN/s data','3000 pN/s data','30000 pN/s data')
xlabel('Bins')
ylabel('Ruptures')

%% 4b. What are the parameters and initial conditions, and are they physiologically possible?

% k = 4.5661e-05;
% fs = 11.1063;
% S1_0 = 46.5254;
% S2_0 = 48.4624;
% S3_0 = 39.5161;
% J = 605.1140

%% 4c. Determine and plot the residuals and/or weighted residuals

resi = data(1:end-1,2:4) - R;                                               %Calculate residuals and plot
plot(data(1:end-1,1),resi);
legend('300 pN/s','3000 pN/s','30000 pN/s')
xlabel('Bins')
ylabel('Ruptures')
title('Residuals')

%% 5a. Address the negative control data

guesses = estimates;                                                        %Set guesses to estimates

[estimates, J] = fminsearch(@objrupture, guesses);                          %Minimize objective function
disp(estimates);

S0 = [estimates(3), estimates(4), estimates(5) estimates(6)];

[f,s] = ode23s(@negODE, fspan, S0, options, estimates, r);                  %Evaluate ODE with estimates

% k = 1.2291e-05
% fs = 10.2486;
% S1_0 = 36.6416;
% S2_0 = 39.7678;
% S3_0 = 31.8495;
% S4_0 = 59.6942;
% a = 0.0060;
% J = 93.3044

R_neg = -diff(s);                                                           %Compute R and plot

R_neg(:,1:3) = R_neg(:,1:3) + R_neg(:,4);

plot(f(1:end-1,1),R_neg,data(:,1),data(:,2),'o',data(:,1),data(:,3),'*',data(:,1),data(:,4),'.',data(:,1),data(:,5),'s');
legend('300 pN/s sim','3000 pN/s sim','30000 pN/s sim','Negative control sim','300 pN/s data','3000 pN/s data','30000 pN/s data','Negative control data')
xlabel('Bins')
ylabel('Ruptures')

close all

resi_neg = data(1:end-1,2:5) - R_neg;                                       %Compute residuals and plot
plot(data(1:end-1,1),resi_neg)
xlabel('Bins')
ylabel('Ruptures')
title('Residuals')

%% 6a. Plot the model and data using a new weighting scheme

guesses = estimates;                                                        %Guess with previous estimates

[estimates, J] = fminsearch(@objrupture, guesses);                          %Minimize objective function
disp(estimates);

S0 = [estimates(3), estimates(4), estimates(5) estimates(6)];

[f,s] = ode23s(@negODE, fspan, S0, options, estimates, r);                  %Evaluate ODE with estimates

R_weight = -diff(s);                                                        %Calculate R and plot
R_weight(:,1:3) = R_weight(:,1:3) + R_weight(:,4);
plot(f(1:end-1,1),R_weight,data(:,1),data(:,2),'o',data(:,1),data(:,3),'*',data(:,1),data(:,4),'.',data(:,1),data(:,5),'s');
legend('300 pN/s sim','3000 pN/s sim','30000 pN/s sim','Negative control sim','300 pN/s data','3000 pN/s data','30000 pN/s data','Negative control data')
xlabel('Bins')
ylabel('Ruptures')

%% 6b. What are the parameters and initial conditions? How are they different from before?

% k = 1.1965e-05;
% fs = 10.2448;
% S1_0 = 39.5175;
% S2_0 = 41.7392;
% S3_0 = 34.0280;
% S4_0 = 45.9888;
% a = 0.0080;
% J = 11.0793

%% 6c. Plot the residuals or weighted residuals. How did the distribution of the residuals or weighted residuals change?

resi_weight = (data(1:end-1,2:5) - R_weight);                               %Calculate residuals and plot
plot(data(1:end-1,1),resi_weight);
xlabel('Bins')
ylabel('Ruptures')
title('Residuals')
