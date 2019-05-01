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

tspan = [0 1];
t_eff = zeros(length(IC_set),4); % num ICs x 4 solvers

for i = 1:length(IC_set)
    fprintf('Working on loop %d/5\n',i)
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
    [t,y] = DSDEexact(@TCellRXN,tspan,IC,params);
    t_eff(i,4) = toc;
end
disp('Finished')
% save("t_eff.mat", 't_eff') % did this once so save it just in case

% plot log time vs. A1_0 and A2_0 value
t_max = max(max(t_eff));

figure;
plot(IC_set,log(t_eff),'-o')
legend('tau @ 0.1','tau @ 0.01','tau @ 0.001','exact')
xlabel('Initial Value of A1,A2')
ylabel('log(time elapsed)')
title('Measuring efficiency of stochastic solvers')

%% Accuracy
clear all; close all; clc

% set A1, A2 IC to 100, use each solver and run as many times in 5 mins
% histogram of A1_p with 10 bins
% fraction of simulations: all on, all on, partially on, error
% avg, std of A1_p at end

% parameters
load("params.mat")
IC(1) = 100;
IC(4) = 100;
tspan = [0 20];

A1_p_end = cell(4,1); % save A1_p end values
A1_p_end{1} = zeros(20000,1);  % cannot be more than this value, determ.
A1_p_end{2} = zeros(10000,1);   % using shortest t_eff for A1 = 100
A1_p_end{3} = zeros(5000,1);    % 1- exact, 2- tau0.1, 3- tau0.01
A1_p_end{4} = zeros(20,1);      % 4- tau0.001
end_index = zeros(4,1);
end_time = 5; % seconds

% for each solver
for solver = 1:4
    i = 0;
    switch solver
        case 1 % exact
            disp('Using exact solver...')
            tic
            while toc < end_time % seconds
                % update index
                i = i+1; 

                % solve for system
                [t,y] = DSDEexact(@TCellRXN,tspan,IC,params);

                % save A1_p(end)
                A1_p_end{solver}(i) = y(end,2);
            end
        case 2 % tau 0.1]
            disp('Using tau 0.1 solver...')
            tic
            while toc < end_time % seconds
                % update index
                i = i+1; 

                % solve for system
                [t,y] = DSDEtauleap(@TCellRXN,tspan,IC,0.1,params);

                % save A1_p(end)
                A1_p_end{solver}(i) = y(end,2);
            end
        case 3 % tau 0.01
            disp('Using tau 0.01 solver...')
            tic
            while toc < end_time % seconds
                % update index
                i = i+1; 

                % solve for system
                [t,y] = DSDEtauleap(@TCellRXN,tspan,IC,0.01,params);

                % save A1_p(end)
                A1_p_end{solver}(i) = y(end,2);
            end
        case 4 % tau 0.001
            disp('Using tau 0.001 solver...')
            tic
            while toc < end_time % seconds
                % update index
                i = i+1; 

                % solve for system
                [t,y] = DSDEtauleap(@TCellRXN,tspan,IC,0.001,params);

                % save A1_p(end)
                A1_p_end{solver}(i) = y(end,2);
            end
    end
end
disp('Finished')

% 
% solver = 1;
% i = 0;
% tic
% while toc < 5 % seconds
%     % update index
%     i = i+1; 
%     
%     % solve for system
%     [t,y] = DSDEexact(@TCellRXN,tspan,IC,params);
%     
%     % save A1_p(end)
%     A1_p_end{solver}(i) = y(end,2);
% end
% % save idx so we know where to stop indexing
% end_index(1) = i;
% disp('Finished')





