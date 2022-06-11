% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Mastrofini Alessandro
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Medical Engineering - University of Rome Tor Vergata
% Physiological Systems Modeling and Simulation
% F. Caselli, MSSF A.Y. 2021/2022
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Glucose minimal model
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Load data
clear all
close all
load("experimental_data.mat")
% data tgi
% time - glucose - insuline
time=tgi(:,1);
glucose=tgi(:,2);
insuline=tgi(:,3);
% plot figure Time [min] vs glucose [mg/dl]
figure;
plot(time,glucose,'ob');
xlabel('Time [min]')
ylabel('Glucose [mg/dl]')
title('Measured plasma glucose concentration')
% plot figure Time [min] vs insuline [uU/ml]
% insulin is measured in Unit and 1 μU/ml = 7.174 pmol/l
figure;
plot(time,insuline,'or');
xlabel('Time [min]')
title('Measured plasma insuline concentration')
ylabel('Insuline [{\mu}U/ml]')

%% parameters for glucose minimal model

G0=279; % [mg/dl]
x0=0;
Gb=93; % baseline glucose concentration
Ib=11;
Sg=2.6E-2;
k=0.025;
Si=5.0e-4; % insuline sensititivty
% color for plotting

mygreen='#77AC30';
myred='#A2142F';
myblue='#0072BD';
%% ODE solver
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

% plot glucose levels
figure;
plot(t,y(:,1),'-','Color',mygreen)
hold on
plot(time,glucose,'o','Color',myred)
legend({'Solved', 'Measured'})
title('Glucose')
xlabel('Time[min]')
ylabel('Glucose [mg/dl]')
legend({'Model','Samples'})

% % plot effective insuline
% figure;
% plot(t,y(:,2),'Color',myblue)
% title('Effective insuline')
% xlabel('Time[min]')
% ylabel('[{\mu}U/ml]')

% plot plasmatic insuline
figure;
time_for_interp=linspace(time(1),time(end),length(time));
I_inter=interp1(2*time,insuline,time_for_interp);
plot(time,I_inter,'-*','Color',myblue)
hold on
plot(time,insuline,'or')
title('Plasmatic levels of insuline')
legend({'Interpolation','Samples'})
xlabel('Time[min]')
ylabel('Insuline [{\mu}U/ml]')
%%
path='figs/';
exportgraphics(figure(1),strcat(path,'glucose_sample','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(2),strcat(path,'insuline_sample','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(3),strcat(path,'glucose_model','.pdf'),'BackgroundColor','none','ContentType','vector');
%exportgraphics(figure(1),strcat(path,'effective_insuline_model','.pdf'),'BackgroundColor','none','ContentType','vector');
exportgraphics(figure(4),strcat(path,'insulin_model','.pdf'),'BackgroundColor','none','ContentType','vector');

%% evaluation
% evaluate solution in the same istant of the samples (time)
sol = ode45(@(t,y) odefcn(t,y,insuline,time,parameters), [time(1), time(end)],[G0,x0]);
% calcolate value G(time), X(time)
evaluated_sol=deval(sol,time);
% calcolate error like (sample-model)./sample
error=100*abs((glucose-evaluated_sol(1,:)')./glucose);
% extract mean error with the first component
error_value=mean(error(2:end));
disp(['Mean error: ',num2str(error_value),' %'])