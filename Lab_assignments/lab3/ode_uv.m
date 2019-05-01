function dFdt = ode_uv(t,y)
    a1 = 3.5;
    a2 = 4;
    
    dFdt = zeros(2,1);
    dFdt(1) = a1 / (1 + y(2)^2) - y(1);
    dFdt(2) = a2 / (1 + y(1)^2) - y(2);
end
