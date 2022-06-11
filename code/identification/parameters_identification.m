% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Mastrofini Alessandro
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Medical Engineering - University of Rome Tor Vergata
% Physiological Systems Modeling and Simulation
% F. Caselli, MSSF A.Y. 2021/2022
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Glucose minimal model
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all
close all
load("experimental_data.mat")
% time=tgi(:,1);
% glucose=tgi(:,2);
% insuline=tgi(:,3);
time=tgi(2:end,1);
glucose=tgi(2:end,2);
insuline=tgi(2:end,3);
%% identification
% parameters
x0=0;
Gb=93; % baseline glucose concentration
Ib=11;
G0=279; % [mg/dl]
Sg=2.6E-2;
Si=5.0e-4;
k=0.025;
% true parameters
p_true=[Sg Si k G0];
% initial guess on parameters
% p_init=[1.5, 1.5, 1.5,1.2].*p_true;
% p_init=p_true;
% from Pacini-Bergman doi: 10.1016/0169-2607(86)90106-9
% p_init=[0.399e-1,0.4e-4,0.2e-1,0.287e3]; 
% p_init=[p1 p2 p3 p4];
% where p1==Sg
% p2==k
% p3==k*S_i
% p4==G0
p1=0.399e-1;
p2=0.2e-1;
p3=0.4e-4;
p4=0.287e3;
p_init=[p1, p3/p2, p2, p4]; 
options = optimoptions('lsqnonlin','Display','iter');
options.Algorithm = 'levenberg-marquardt';
options.StepTolerance=1e-8;
options.PlotFcn='optimplotresnorm';
options.Display='iter-detailed';
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@optfcn,p_init,[],[],options);

disp(['Error: ',num2str((p_true-x)./p_true)])
disp(['p_true: ' , num2str(p_true)])
disp(['p_optm: ' , num2str(x)])

% plot results vs sample data
G0=x(4);
x0=0;
Gb=93; % baseline glucose concentration
Ib=11;
Sg=x(1);
Si=x(2);
k=x(3);
parameters=[Sg,Gb,k,Ib,Si]; 
[t,y] = ode45(@(t,y) odefcn(t,y,insuline,time,parameters), time,[G0,x0]);
%figure;plot(t,y(:,1));hold on;plot(time,glucose,'-o')
figure;plot(t,y(:,1));hold on;plot(time,glucose,'-o')
legend({'Identification','Samples'})
xlabel('Time [s]')
ylabel('Glucose concentration [mg/dl]')

%%
exportgraphics(figure(2),'optim.pdf','BackgroundColor','none','ContentType','vector');