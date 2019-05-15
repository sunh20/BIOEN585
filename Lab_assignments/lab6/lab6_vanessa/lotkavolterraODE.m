function dydt = lotkavolterraODE(t,y,p)
% ode for the model:
dydt = [ p(1)*y(1)-p(2)*y(1)*y(2)         
         p(3)*y(1)*y(2)-p(4)*y(2) ];
