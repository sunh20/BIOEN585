%% Lab 1
% Samantha Sun
% BIOEN 585
% 20180404

clear all; close all; clc

%% Question 1: Electrical Model
% first-order ODE, solve numerically for I
clearvars

% define constants
r1 = 0.2; % ohm
r2 = 0.2; % ohm
c = 1; % farad

% time step + time values
dt = 0.01; % seconds
t = 0:dt:2; 

% variables
I = zeros(1,length(t));
I(1) = 0; % intial condition
V = zeros(1,length(t)); 
V = V + 5; % for t > 0, V = 5, does not vary with time
V(1) = 0; 

% step through loop
for i = 1:length(t)-1
    I(i+1) = I(i) + dt * (1/c/r1/r2) * (V(i) - I(i)*(r1+r2));
end

% plot response
figure;
plot(t, I)
xlabel('Time (s)')
ylabel('Current (Amp)')
title('Response curve of current (I) in model')

% change voltage generator input
I2 = zeros(1,length(t));
V2 = 5*sin(10*t); % given from assignment

% step through loop
for i = 1:length(t)-1
    I2(i+1) = I2(i) + dt * (1/c/r1/r2) * (V2(i) - I2(i)*(r1+r2));
end

% plot response
figure;
subplot(2,1,1)
plot(t, I2)
xlabel('Time (s)')
ylabel('Current (Amp)')
title('Response curve of current (I) in model')

subplot(2,1,2)
plot(t, V2)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Input voltage function')

%% Question 2: Fluid model
% solve for flow in flow chamber

clearvars

% time
dt = 0.01;
t = -1:dt:7;

% parameters
Rf = 5e9;   % Pa*s/m^3
Ct = 1e-10; % m^3/Pa
Q = zeros(1, length(t));
Q(t > 0 & t < 3) = 2e-9; % flow starts at t=0, stops at t=3

% variables
Qf = zeros(1,length(t));

% numeric solution
for i = 1:length(t)-1
    Qf(i+1) = Qf(i) + dt * ((Q(i)-Qf(i))/Rf/Ct);
end

% plot
figure;
plot(t, Q, t, Qf)
xlabel('Time (s)')
ylabel('Flow Rate (m^3/s)')
title('Response curve of flow rates (Q) in model')
legend('Input Pump','Flow Chamber')
ylim([0 2.2e-9])

%% Question 2: Fluid model using ode function
% solve for flow in flow chamber

clearvars

% time
dt = 0.01;
t = -1:dt:7;

% parameters
Rf = 5e9;   % Pa*s/m^3
Ct = 1e-10; % m^3/Pa
Q = zeros(length(t),1);
Q(t > 0 & t < 3) = 2e-9; % flow starts at t=0, stops at t=3

% variables
Qf = zeros(length(t),1);

% numeric solution
for i = 1:length(t)-1
    Qf(i+1) = Qf(i) + dt * ((Q(i)-Qf(i))/Rf/Ct);
end

% use ode45 function
[t2, Qf2] = ode45(@(t,y) q2ode(t,y,Rf,Ct), t, 0);

% use ode stiff solver function
[t3, Qf3] = ode15s(@(t,y) q2ode(t,y,Rf,Ct), t, 0);

% change abs tolerance
options = odeset('AbsTol',1e-15);
[t4, Qf4] = ode45(@(t,y) q2ode(t,y,Rf,Ct), t, 0, options);

% plot
figure;
plot(t, Qf, t2, Qf2, t3, Qf3, t4, Qf4)
xlabel('Time (s)')
ylabel('Flow Rate (m^3/s)')
title('Comparing numerical vs ode function solution')
legend('Numerical','ode45','ode15s','ode45+change abs tol')
ylim([0 2.2e-9])


%% Question 3: Chemical Model
% solve for concentrations of A1, A2 over time
% M = moles

% time
dt = 10;
t = 0:dt:1000;

% variables
A1 = zeros(1,length(t));
A2 = zeros(1,length(t));
C2 = zeros(1,length(t));

% parameters
D = 20e-6;      % M
V1 = 1;         % L
V2 = 5;         % L
RT = 5e-6;      % M
k_off = 0.001;  % s^-1
k_on = 1e4;     % s^-1 M^-1
P1 = 0.01;      % L/s
P2 = 0.0005;    % L/s

% initial conditions
A1(1) = D;

% solve
for i = 1:length(t)-1
    A1(i+1) = A1(i) + dt * (-P1*A1(i)/V1);
    A2(i+1) = A2(i) + dt * (P1*A1(i)/V2 - P2*A2(i)/V2 + k_off*C2(i) - ...
                                                k_on*A2(i)*(RT-C2(i)));
    C2(i+1) = C2(i) + dt * (k_on*A2(i)*(RT-C2(i)) - k_off*C2(i));
end

% plot
figure;
subplot(3,1,1)
plot(t,A1*1e6)
xlabel('Time (s)')
ylabel('Concentration (umol/L)')
title('Response curve of modeled drug concentration in stomach')

subplot(3,1,2)
plot(t,A2*1e6)
xlabel('Time (s)')
ylabel('Concentration (umol/L)')
title('Response curve of modeled drug concentration in blood')

subplot(3,1,3)
plot(t,C2*1e6)
xlabel('Time (s)')
ylabel('Concentration (umol/L)')
title('Response curve of modeled receptor-drug complex concentration')

%% Question 4: Mechanical model
clearvars

% time
dt = 0.01;
t = 0:dt:10;

% variables
x = zeros(length(t),1);

% params
b = 0.8;    % Ns/m
k = 1;      % N/m
m = 0;      % kg
F = 3e-3;   % N

% initial condition, x = 0
% starting muscle length = 2.5mm, add on to end
% because we're solving for change in mucle length

% solve
for i = 1:length(t)-1
    x(i+1) = x(i) + dt * (F - k*x(i)) / b;
end

% make adjustments
x = x * 1e3; % convert to mm
x = x + 2.5; % muscle length @ equilibrium

% plot
figure;
plot(t, x)
ylim([0 6])
xlabel('Time (s)')
ylabel('Length (mm)')
title('Modeled change in muscle length')












