% Calculate the big term from rho_0 equation.
function bterm=bigterm(bet,lam,t_L,t_C)
    bterm = 0;
    for i = 1:6
        bterm = bterm + bet(i)/(1.0 + ((1.0-exp(-lam(i)*t_L))/(lam(i)*t_C)));
    end
end