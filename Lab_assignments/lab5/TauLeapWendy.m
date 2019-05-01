function [t,X]=TauLeapWendy(DefineReactions,tspan, IC, RelTol, params)
%Tau Leap Algorithm for chemical reaction equations
%BY WENDY THOMAS
% S = stoichiometry of C substrates in R reactions
% P = stoichiometry of products
% K = vector of reaction rates

% INITIALIZE by calling the passed function to define the reactions:
if ~isa(DefineReactions, 'function_handle') ...
        disp('pass a function handle with your reaction defs'); return; end
X = IC;  % initialize the matrix of chemicals at time t
[S, P, K] = DefineReactions(params); % defines reaction stoichiometry
% check that size of inputs are correct:
[R, C] = size(S);
if size(S) ~= size(P); disp('reaction defs inconsistent'); return;end
if length(IC) ~= C; disp('ICs inconsistent with reactions'); return; end
if size(K,2) == 1; K = K'; end  % force column vector
if length(K) ~= R; disp ('parameters inconsistent with reactions');return;end
if length(tspan) ~= 2; disp('Time vector: [0;maxtime] needed'); return;end;
maxT = tspan(2);
tau = maxT/1000; % initial check; time step is dynamic.
rand('state',sum(100*clock));  %set the random number generator
t = 0;  %time
% run loop until time is up or all chemicals are gone:
while (t(end)<maxT && any(X(end,:)) )
    %step 1: Calculate a_r and lambda(r)
    a = K;
    for r = 1:R;  
        for c = 1:C
            if S(r,c) == 1; a(r) = a(r)*X(end,c);
            elseif S(r,c) == 2; 
                a(r) = a(r)*X(end,c)*(X(end, c)-1)/2;
            elseif S(r,c) == 3; 
                a(r) = a(r)*X(end,c)*(X(end, c)-1)/2*(X(end,c)-2)/3;
            end; 
        end
    end 
    if a == 0 ;% system can't change; jump to end and exit;
        X(end+1,:) = X(end,:);
        t =[t; maxT];
        break; 
    end 
    %Step 2a; test time step for right size
    % do this by testing the maximum effect of the expected reactions:    
    tau = tau*2; % always try to increase time step so not stuck short;
    lambda = a*tau;
    Xcurrent = max(X(end,:),ones(1,C));
    Max_change = max(abs(lambda*(-S + P)./Xcurrent));
    % but then shorten time step if needed until within tolerance:
    while Max_change > RelTol;
        tau = tau/2;
        lambda = a*tau;
        Max_change = max(abs(lambda*(-S + P)./Xcurrent));
    end
    %Step 2: determine how much of each reaction occurs:
    lambda = a*tau; % expected number
    for r = 1:R; 
        if lambda(r) < 20 
            k(r) = random('poisson',lambda(r));
        else k(r) = round(lambda(r)+sqrt(lambda(r))*randn);
            if k(r) < 0; k(r) = 0; end 
            % make sure normal approximation never makes negative rxns
        end
    end
    %Step 3: carry out the reaction r
    t =[t; t(end) + tau]; % t is time array; add last entry to it.
    new = X(end,:)+ k*(-S + P);
    X = [X ; new];
end