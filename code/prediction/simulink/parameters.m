G0=279; % [mg/dl]
x0=0;
Gb=93; % baseline glucose concentration
Ib=11;
Sg=2.6E-2;
k=0.025;
Si=5.0e-4;
load("experimental_data.mat")
% data tgi
% time - glucose - insuline
% plot figure Time [min] vs glucose
time=tgi(:,1);
glucose=tgi(:,2);
insuline=tgi(:,3);