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

% plotting polynomial
u = 0:0.01:5;
v = 0:0.01:5;
u_null = a1./(1+v.^2);
v_null = a2./(1+u.^2);

figure;
plot(u,v_null)
hold on
plot(u_null,v)
plot(null_u, null_v, 'o')
xlabel('u')
ylabel('v')
xlim([0 5])
ylim([0 5])
title('Phase Plane Plot of u, v, nullclines + equilibrium points')

% stability analysis - get eigenvalues @ each eq point
J1 = [-1, -2*null_v(1)*a1/((1+null_v(1)^2)^2);
        -2*null_u(1)*a2/((1+null_u(1)^2)^2), -1];
[V1,D1] = eig(J1);
    
J2 = [-1, -2*null_v(2)*a1/(1+null_v(2)^2)^2;
    -2*null_u(2)*a2/(1+null_u(2)^2)^2, -1];
[V2,D2] = eig(J2);

J3 = [-1, -2*null_v(3)*a1/(1+null_v(3)^2)^2;
    -2*null_u(3)*a2/(1+null_u(3)^2)^2, -1];
[V3,D3] = eig(J3);

%% Question 2: Phase Plane Analysis

clearvars -except null_u null_v a1 a2 V1 D1 V2 D2 V3 D3

% define eq points
eq1 = [null_u(1), null_v(1)];
eq2 = [null_u(2), null_v(2)];
eq3 = [null_u(3), null_v(3)];

% solve trajectories for array of initial conditions

% real-time plotting magic
figure;
plot(null_u, null_v, 'o')
hold on
xlabel('u')
ylabel('v')
xlim([0 5])
ylim([0 5])
title('Phase Plane Plot of u, v trajectories')

tspan = [0 50];
for i = 0:5
    u0 = i;
    for j = 0:5
        v0 = j;
        [t,y] = ode45(@ode_uv,tspan,[u0 v0]);
        u = y(:,1);
        v = y(:,2);
        
        % colors + add to plot
        if [round(u(end),2), round(v(end),2)] == round(eq1,2)
            plot(u,v,'Color','[0.929 0.694 0.098]') % gold
        end
        
        if [round(u(end),2),round(v(end),2)] == round(eq2,2)
            plot(u,v,'c') % cyan
        end
        
        if [round(u(end),2),round(v(end),2)] == round(eq3,2)
            plot(u,v,'Color','[0.494 0.184 0.556]') % PURPLE
        end
        %pause(0.1)
        
    end
end

% plotting eigenvalues + eigenvectors (EXTRA CREDIT)
% line centered @ each EQ point, length of eigenvalue + direction
% of eigenvector 
% color black or red depending on sign of eigenvalue
% black = pos, red = neg
Ds = [diag(D1); diag(D2); diag(D3)];
EQs = [eq1; eq2; eq3];
Vs = [V1; V2; V3];

for idx = 1:3
    D = Ds(2*idx-1:2*idx);
    eq = EQs(idx,:);
    V = Vs(2*idx-1:2*idx,:);
    x_eq1 = [eq(1) - V(1,1)*sqrt(2)*D(1),eq(1) + V(1,1)*sqrt(2)*D(1)];
    y_eq1 = [eq(2) - V(2,1)*sqrt(2)*D(1),eq(2) + V(2,1)*sqrt(2)*D(1)];
    if D(1) > 0
        plot(x_eq1,y_eq1,'-k')
    else
        plot(x_eq1,y_eq1,'-r')
    end

    x_eq1b = [eq(1) - V(1,2)*sqrt(2)*D(2),eq(1) + V(1,2)*sqrt(2)*D(2)];
    y_eq1b = [eq(2) - V(2,2)*sqrt(2)*D(2),eq(2) + V(2,2)*sqrt(2)*D(2)];
    plot(x_eq1b,y_eq1b,'-r')
    
    if D(1) > 0
        plot(x_eq1b,y_eq1b,'-k')
    else
        plot(x_eq1b,y_eq1b,'-r')
    end
end

%% Question 3: Bifurcation Analysis
