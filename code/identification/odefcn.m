function dydt=odefcn(t,y,insuline,time,parameters)
% exctract parameters from the vector

Sg=parameters(1);
Gb=parameters(2);
k=parameters(3);
Ib=parameters(4);
Si=parameters(5);

% define ODE function to solve

% first I need to reconstruct I(t) 
% the insuline plasmatic concentration
% ODE solver need I(t) at every t
% so with interp1() is possibile to interpolate
% the insuline (from file) value to the correspondent 
% value of t
% interpolate in the same variabile of the ode solver (t)
I_inter=interp1(time,insuline,t);

% now it's possibile to define function 
% y(1) = G(t)
% y(2)= X(t)
% system of the ODE
% dGdt=Sg (Gb-G(t))-X(t)*G(t)
% dXdt=k*(Si*(I(t)-Ib)-X(t)) 
% NOTE: I(t) isn't insuline variable and must be interpolated
dydt(1)=Sg*(Gb-y(1))-y(2)*y(1);
dydt(2)=k*(Si*(I_inter-Ib)-y(2));
dydt=dydt'; % trasponse to ODE45 sintax
end


