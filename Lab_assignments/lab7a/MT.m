function [X,Y,Q] = MT(params,tspan,dt,N,d)
% solves for N number of microtubule tracks over time t
% wall diameter d if wall is present

if ~exist('d','var')
    d = NaN;
end

% define needed terms
Dv = params(1);
V_avg = params(2);
Lp = params(3);

R_std = sqrt(2*Dv*dt);
Q_std = sqrt(V_avg*dt/Lp);

% variables - IC = 0
X = nan(length(tspan),N);
Y = nan(length(tspan),N);
Q = nan(length(tspan),N);

X(1,:) = 0;
Y(1,:) = 0;
Q(1,:) = 0;

% calculate tracks
for t = 2:length(tspan)
    dQ = normrnd(0,Q_std,[1,N]);
    r = normrnd(V_avg*dt,R_std,[1,N]);
    
    Q(t,:) = Q(t-1,:) + dQ;     % new sigma term
    dx = r.*cos(Q(t,:));
    dy = r.*sin(Q(t,:));
    
    X(t,:) = X(t-1,:) + dx;
    Y(t,:) = Y(t-1,:) + dy;
    
    % recalculate if there is a wall
    idxs = abs(Y(t,:)) > (d/2); % indices where particle passes bounds
    
    if sum(idxs) > 0
        disp('I hit a wall')
        Q(t,:) = 0;                     % reset to zero
        %X(t,idxs) = X(t-1,idxs) + (d/2-Y(t-1,idxs)) ./ ...
        %        sin(Q(t-1,idxs)) .* (cos(Q(t-1,idxs))-1) + r(idxs); % trig
        X(t,idxs) = X(t-1,idxs) + r(idxs);
        Y(t,idxs) = sign(Y(t-1,idxs)).*d/2; % cling to border

        % check that it was done correctly
        if sum(abs(Y(t,:)) > (d/2)) > 1
            disp('something went wrong')
        end  
    end
    
%     % plot to watch tracks grow!
%     plot(X(1:t,:),Y(1:t,:)) 
%     if ~isnan(d) == 1
%         hold on
%         plot(tspan(1:t),ones(1,t)*d/2,'k--','LineWidth',2)
%         plot(tspan(1:t),ones(1,t)*-d/2,'k--','LineWidth',2)
%     end
%     xlabel('x position')
%     ylabel('y position')
%     title('Microtubule motion')
%     drawnow
%     pause(0.01)
end


end