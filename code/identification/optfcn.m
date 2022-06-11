function X=optfcn(x)

load("experimental_data.mat");
% time=tgi(:,1);
% glucose=tgi(:,2);
% insuline=tgi(:,3);

time=tgi(2:end,1);
glucose=tgi(2:end,2);
insuline=tgi(2:end,3);

% [G0 Sg Si k]
G0=x(4);
x0=0;
Gb=93; % baseline glucose concentration
Ib=11;
Sg=x(1);
Si=x(2);
k=x(3);
% ODE solver
parameters=[Sg,Gb,k,Ib,Si]; 

%[t,y] = ode45(@(t,y) odefcn(t,y,insuline,time,parameters), time,[G0,x0]);

sol = ode45(@(t,y) odefcn(t,y,insuline,time,parameters), [time(1), time(end)],[G0,x0]);
% calcolate value G(time), X(time)
evaluated_sol=deval(sol,time);
% error function
X=glucose-evaluated_sol(1,:)';
end