% Lezione n. 10 05/04/2022

% Load data from .mat file
load("experimental_data.mat")
% data tgi
% time - glucose - insuline
% plot figure Time [min] vs glucose
time=tgi(:,1);
glucose=tgi(:,2);
insuline=tgi(:,3);

figure;
plot(time,glucose,'ob');
xlabel('Time [min]')
ylabel('[mg/dl]')
title('Measured plasma glucose concentration')
figure;
plot(time,insuline,'or');
xlabel('Time [min]')
title('Measured plasma insuline concentration')
ylabel('[mU/ml ???]')

%%
% abbiamo estratto le quantità I(t) 
% calcolare G(t) e X(t)

% parametri:

G0=279; % [mg/dl]
x0=0;
Gb=93; % baseline glucose concentration
Ib=11;
Sg=2.6E-2;
k=0.025;
Si=5.0e-4; % insuline sensititivty
%%
parameters=[Sg,Gb,k,Ib,Si]; % insert all the parameters inside a vector
% call ODE45 by passing odefnc() 
% passe odefnc() with t,y and insuline (must be interpolated)
% and time (from tgi variable) and parameters vector
% use span time as a [init final] and not all the time vector
% so ODE45 estimates y in t value by itself 
% and not in the same value of the time variabile
% So, using time variabile ODE return y with 24 components
% using [time(1) time(end)] ODE45 return y in 69 components

[t,y] = ode45(@(t,y) odefcn(t,y,insuline,time,parameters), [time(1), time(end)],[G0,x0]);

figure;
plot(t,y(:,1),'-r')
hold on
plot(time,glucose,'ob')
legend({'Solved', 'Measured'})
title('Glucose')
xlabel('Time[min]')
ylabel('[mg/dl]')
figure;
plot(t,y(:,2))
title('Effective insuline')


