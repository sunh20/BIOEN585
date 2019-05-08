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

% % Stochastic soln: exact
% tspan = [0 20];
% [t,y] = DSDEexact(@TCellRXN,tspan,IC,params);

% Stochastic soln: tau 0.1
tspan = [0 20];
[t,y] = DSDEtauleap(@TCellRXN,tspan,IC,0.1,params);
y(end,2)

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
A1_p_end{1} = zeros(15000,1);  % cannot be more than this value, determ.
A1_p_end{2} = zeros(10000,1);   % using shortest t_eff for A1 = 100
A1_p_end{3} = zeros(5000,1);    % 1- exact, 2- tau0.1, 3- tau0.01
A1_p_end{4} = zeros(20,1);      % 4- tau0.001
end_index = zeros(4,1);
end_time = 5; % seconds

% for each solver
for solver = 1:3 % skip tau0.001 because it takes too long
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
                A1_p_end{solver}(i) = y(end-2,2);
            end
            % save idx so we know where to stop indexing
            end_index(solver) = i;
        case 2 % tau 0.1]
            disp('Using tau 0.1 solver...')
            tic
            while toc < end_time % seconds
                % update index
                i = i+1; 

                % solve for system
                [t,y] = DSDEtauleap(@TCellRXN,tspan,IC,0.1,params);

                % save A1_p(end)
                A1_p_end{solver}(i) = y(end-2,2);
            end
            % save idx so we know where to stop indexing
            end_index(solver) = i;
        case 3 % tau 0.01
            disp('Using tau 0.01 solver...')
            tic
            while toc < end_time % seconds
                % update index
                i = i+1; 

                % solve for system
                [t,y] = DSDEtauleap(@TCellRXN,tspan,IC,0.01,params);

                % save A1_p(end)
                A1_p_end{solver}(i) = y(end-2,2);
            end
            % save idx so we know where to stop indexing
            end_index(solver) = i;
        case 4 % tau 0.001
            disp('Using tau 0.001 solver...')
            tic
            while toc < end_time % seconds
                % update index
                i = i+1; 

                % solve for system
                [t,y] = DSDEtauleap(@TCellRXN,tspan,IC,0.001,params);

                % save A1_p(end)
                A1_p_end{solver}(i) = y(end-2,2);
            end
            % save idx so we know where to stop indexing
            end_index(solver) = i;
    end
end
disp('Finished')
% save("A1_p_end.mat",'A1_p_end','end_index')
%% post processing
load("A1_p_end.mat")
% trim each cell
for i = 1:3
    A1_p_end{i} = A1_p_end{i}(1:end_index(i));
end

% plot histogram of A1_prot w/ 10 bins
str = ["exact";"tau0.1";"tau0.01"];
figure;
for i = 1:3 
    subplot(3,1,i)
    histogram(A1_p_end{i},10) % 10 bins
    xlim([0 100])
    title(sprintf('End value of A1-prot using %s solver',str(i)))
end

% calculate fraction of simulation results + stats
% 1- all on, 2- all off, 3- partially on, 4- error (NaN, negative...etc)
% 5- avg, 6- std
%rows = results, cols = simulators
results = zeros(4,6);
for i = 1:3
    dat = A1_p_end{i};
    results(i,1) = length(dat(dat == 100)); 
    results(i,2) = length(dat(dat == 0)); 
    results(i,3) = length(dat(dat > 0 & dat < 100)); 
    results(i,4) = length(dat) - results(i,1) - results(i,2) - results(i,3); 
    results(i,5) = mean(dat);
    results(i,6) = std(dat);
end

%% make pretty table
AllOn = results(1:3,1);
AllOff = results(1:3,2);
PartOn = results(1:3,3);
Error = results(1:3,4);
Mean = results(1:3,5);
StdDev = results(1:3,6);

table(str,AllOn,AllOff,PartOn,Error,Mean,StdDev)

%% Q3: Run @ different ICs to see how system changes
% use exact solver
clear all; close all; clc

% parameters
load("params.mat")
IC_set = [10, 100, 1000, 10000, 100000]; % for A1_0, A2_0
tspan = [0 20];

% save A1_p_end for each condition
A1_p_end2 = cell(length(IC_set),1);
A1_p_end2{1} = zeros(100000,1);
A1_p_end2{2} = zeros(5000,1);
A1_p_end2{3} = zeros(1000,1);
A1_p_end2{4} = zeros(1000,1);
A1_p_end2{5} = zeros(1000,1);

end_time = 300; % 5 mins
end_index = zeros(1,length(IC_set));

for i = 1:length(IC_set)
    fprintf('Working on loop %d/5\n',i)
    IC(1) = IC_set(i);
    IC(4) = IC_set(i);
    
    ndIdx = 0;

    tic
    while toc < end_time % seconds
        % update index
        ndIdx = ndIdx+1; 

        % solve for system
        if i < 3
            [t,y] = DSDEexact(@TCellRXN,tspan,IC,params);
        else
            [t,y] = DSDEtauleap(@TCellRXN,tspan,IC,0.1,params);
        end

        % save A1_p(end)
        A1_p_end2{i}(ndIdx) = y(end-2,2);
    end
    % save idx so we know where to stop indexing
    end_index(i) = ndIdx;
end
disp('Finished')

save("A1_p_end2.mat",'A1_p_end2','end_index')
%% post processing
load("A1_p_end2.mat")

% trim each cell
for i = 1:length(end_index)
    A1_p_end2{i} = A1_p_end2{i}(1:end_index(i));
end

% plot histogram of A1_prot w/ 10 bins
str = ["A1 = 1e1"; "A1 = 1e2"; "A1 = 1e3"; "A1 = 1e4"; "A1 = 1e5"];
figure;
for i = 1:length(end_index)
    subplot(length(end_index),1,i)
    histogram(A1_p_end2{i},10) % 10 bins
    xlim([0 10^i])
    title(sprintf('End value of A1-prot with IC %s',str(i)))
end

% calculate fraction of simulation results + stats
% 1- avg frac prot, 2- std frac prot, 3- all on
%rows = results, cols = simulators
IC_set = [10, 100, 1000, 10000, 100000];
results = zeros(length(end_index),3);
for i = 1:length(end_index)
    dat = A1_p_end2{i};
    dat_frac = dat/IC_set(i);
    results(i,1) = mean(dat_frac);
    results(i,2) = std(dat_frac); 
    results(i,3) = length(dat_frac(dat_frac == 1));  
end

% make pretty table
Mean = results(:,1);
StdDev = results(:,2);
AllOn = results(:,3);

table(str,Mean,StdDev,AllOn)






















