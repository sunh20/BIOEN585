function dydt = Tcellode(t,y,params)

k1 = params(1);
k2 = params(2);
k3 = params(3);
k4 = params(4);
k5 = params(5);
kD = params(6);

A1 = y(1);
A1_p = y(2);
A1_i = y(3);
A2 = y(4);
E = y(5);
S = y(6);

dA1 = - k3*E*A1 - k5*S*A1;
dA1_p = k3*E*A1;
dA1_i = k5*S*A1;
dA2 = 0;
dE = k1*A1 + k4*A1_p - kD*E;
dS = k2*A2 - kD*S;

dydt = [dA1;dA1_p;dA1_i;dA2;dE;dS];

end