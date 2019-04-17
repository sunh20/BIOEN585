%% Lab 2: Linear Systems
% Samantha Sun
% BIOEN 585
% 20190410

clear all; close all; clc
plotOn = 0;

% Modeling Lung Mechanics
% Pressure-input P(t)
% Measure volume of air in lungs V(t)

% variables
Rc = 0.001; % atm*s/L
Rp = 0.0005;
Cl = 200;   % L/atm
Cw = 200;
Cs = 5;
Ct = (1/Cl+1/Cw+1/Cs)^-1;

% transfer function
num = [1, 1/(Rp*Ct)];
den = [Rc, 1/Cs+Rc/Rp/Ct, 1/Rp/Cs*(1/Cl+1/Cw)];
H = tf(num,den);

% periodic input
maxP = 4.76e-3;     % atm
T = 3;              % s
t = 0:0.001:10;     % time vector
input = square(2*pi*t/T+pi);
input = (input+1)/2*maxP; % adjust to params
resp = lsim(H,input,t);

if plotOn == 1
    % plot response to step function
    figure;
    step(H);
    xlabel('Time')
    ylabel('Volume (L)')
    title('System response to step function')
    xlim([0 2])
    
    % plot response to periodic input
    figure;
    plot(t,resp)
    xlabel('Time')
    ylabel('Volume (L)')
    title('System response to periodic function')
end

% calculate tidal volume (difference between min + max)
TV = max(resp) - min(resp);
fprintf('Tidal volume: %0.2d Liters \n',TV)

%% state space model - steady state gain
A = [-200, 200; 400, -420];
B = [1000; 0];
C = 1;
D = 0;

SSG = -C*A^-1*B;

%% Diseased patients - decreased compliance of alveoli
Cl_mod = [200, 100, 40, 20];
resp_mod = [];

for i = 1:4
    Ct = (1/Cl_mod(i)+1/Cw+1/Cs)^-1;
    num = [1, 1/(Rp*Ct)];
    den = [Rc, 1/Cs+Rc/Rp/Ct, 1/Rp/Cs*(1/Cl_mod(i)+1/Cw)];
    H_mod = tf(num,den); 
    resp_mod = [resp_mod, lsim(H_mod,input,t)];
end

figure;
plot(t,resp_mod)
xlabel('Time')
ylabel('Volume (L)')
title('System response to periodic function - Healthy vs. non-healthy')
legend('Healthy','2-fold','5-fold','10-fold')


%% Diseased patients - increase pressure protocol
input_mod = [1 1.5 2.7 4.5];
resp_mod_help = [];

for i = 1:4
    Ct = (1/Cl_mod(i)+1/Cw+1/Cs)^-1;
    num = [1, 1/(Rp*Ct)];
    den = [Rc, 1/Cs+Rc/Rp/Ct, 1/Rp/Cs*(1/Cl_mod(i)+1/Cw)];
    H_mod = tf(num,den); 
    resp_mod_help = [resp_mod_help, lsim(H_mod,input*input_mod(i),t)];
end

figure;
a1 = subplot(2,1,1);
plot(t,resp_mod)
xlabel('Time')
ylabel('Volume (L)')
title('System response to periodic function - Healthy vs. non-healthy')
legend('Healthy','2-fold','5-fold','10-fold')

a2 = subplot(2,1,2);
plot(t,resp_mod_help)
xlabel('Time')
ylabel('Volume (L)')
title('Non-healthy + increased pressure')
legend('x1','x1.5','x2.7','x4.5')

linkaxes([a1 a2],'xy')

%% fractional sensitivity of SSG, slowest time constant, tidal volume
% to compliance of alveoli

Rc = 0.001; % atm*s/L
Rp = 0.0005;
Cw = 200;
Cs = 5;
Cl = 1:200;   % L/atm
Ct = (1./Cl+1/Cw+1/Cs).^-1;

num = zeros(200,2);
num(:,1) = 1;
num(:,2) = 1./(Rp*Ct);

den = zeros(200,3);
den(:,1) = Rc;
den(:,2) = 1/Cs+Rc/Rp./Ct;
den(:,3) = 1/Rp/Cs*(1./Cl+1/Cw);

% periodic input
maxP = 4.76e-3;     % atm
T = 3;              % s
t = 0:0.001:3;     % time vector
input = square(2*pi*t/T+pi);
input = (input+1)/2*maxP; % adjust to params

% steady state gain
SSG = (Cs*Cw + Cs*Cl + Cl*Cw)./(Cl + Cw);

% shortest time constant tidal volume
sT = zeros(200,1);
TV = zeros(200,1);

for i = 1:length(sT)
    sT(i) = -1/max(roots(den(i,:)));
    
    H = tf(num(i,:),den(i,:));
    resp = lsim(H,input,t);
    TV(i) = max(resp) - min(resp);

end

% plots
figure; 
ax1 = subplot(3,1,1);
plot(Cl, SSG)
xlabel('Alveoli Compliance (L/atm)')
ylabel('Steady State Gain (L/atm)')
title('Change in Alveoli Compliance vs. Outcomes')

ax2 = subplot(3,1,2);
plot(Cl, sT)
xlabel('Alveoli Compliance (L/atm)')
ylabel('Shortest time constant (s)')

ax3 = subplot(3,1,3);
plot(Cl, TV)
xlabel('Alveoli Compliance (L/atm)')
ylabel('Tidal Volume (L)')

linkaxes([ax1, ax2, ax3], 'x')
%TV = max(resp) - min(resp);
