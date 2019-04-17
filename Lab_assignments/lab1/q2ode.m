function dydt = q2ode(t,Qf2,Rf,Ct)
	if (t>0 && t<3)
        Q = 2e-9;
    else
        Q = 0;
    end
    dydt = (Q-Qf2)/Rf/Ct;
end