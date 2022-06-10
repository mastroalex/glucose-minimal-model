% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Mastrofini Alessandro
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Medical Engineering - University of Rome Tor Vergata
% Physiological Systems Modeling and Simulation
% F. Caselli, MSSF A.Y. 2021/2022
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Glucose minimal model
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

load("../experimental_data.mat")
time=tgi(:,1);
glucose=tgi(:,2);
insuline=tgi(:,3);

x0=0;
Gb=93; % baseline glucose concentration
Ib=11;

G0=279; % [mg/dl]
Sg=2.6E-2;

Si=5.0e-4;
k=0.025;
% dobbiamo trovare Si, Sg, G0, K
% 
p_true=[G0 Sg Si k];
% initial guess on parameters
p_init=[0.5, 1.5, 0.2, 2].*p_true;
%p_init=p_true;

%%%% TRY MORE AND TRY TO FIND DIVERGENGE:
%p_init=[4; 7].*p_true;

%% 
% relaxation parameter
% alpha=1;
alpha=0.1;

% regularization parameter
% lambda=0;
lambda=0.05;

% tolerance
tol=1e7*eps;

% maximum number of iterations
max_iter=35;

% iterazione
done=false;

% initialize parameters
p=p_init;

% enable/disable parameter convergence plot
plot_flag=true;
max_iter_plot=10;

if plot_flag
    % plot true parameter value and current guess
    f1=figure();
    subplot(2,2,1)
    hold on
    line([0,max_iter_plot],[p_true(1),p_true(1)],'linewidth',2,'Color','r','LineStyle',':')
    xlabel('Iteration')
    ylabel('Parameter G0')
    set(gca,'FontSize',12)

    subplot(2,2,2)
    hold on
    line([0,max_iter_plot],[p_true(2),p_true(2)],'linewidth',2,'Color','b','LineStyle',':')
    xlabel('Iteration')
    ylabel('Parameter Sg')
    set(gca,'FontSize',12)

     subplot(2,2,3)
    hold on
    line([0,max_iter_plot],[p_true(3),p_true(3)],'linewidth',2,'Color','b','LineStyle',':')
    xlabel('Iteration')
    ylabel('Parameter Si')
    set(gca,'FontSize',12)

     subplot(2,2,4)
    hold on
    line([0,max_iter_plot],[p_true(4),p_true(4)],'linewidth',2,'Color','b','LineStyle',':')
    xlabel('Iteration')
    ylabel('Parameter K')
    set(gca,'FontSize',12)
end

i=0;
% Gauss-Newton iteration
while (~done) && (i<=max_iter)
    i=i+1;
     %parameters=[Sg,Gb,k,Ib,Si];
    % p_true=[G0 Sg Si k];
   parameters=[p(2), Gb, p(4), Ib, p(3)];
   [t,y] = ode45(@(t,y) odefcn(t,y,insuline,time,parameters), time,[p(1),x0]);

     y_pred=y(:,1);
    % objective funtion to be minimized: E=1/2*||e||^2
    disp('ODE concluso')
    % error vector
    e=(y_pred-glucose); % *
    
    % Jacobian matrix - sensitivity wrt parameter (updated at each iteration)
    J=jacobian_fun(p,time,Gb,Ib,x0,insuline);
    disp('J fun concluso')

    % gradient of objective function wrt parameters
    dE_dp=J'*e;
    
    % Hessian of objective function with second derivatives neglected (Gauss-Newton approximation)
%     d2E_dp2=J'*J;
    % % regularized Hessian
    d2E_dp2=J'*J+lambda*eye(4);
    
    h=-d2E_dp2\dE_dp;
    
    % check (stop criterion on h)
    done=norm(h)<tol;
    
    p_save=p; % for plotting purpose
    
    % parameter update (with relaxation)
    p=p+alpha*h';
    
    if plot_flag
        % update parameter convergence plot
        figure(f1);
        subplot(2,2,1)
        hold on
        plot([i-1,i],[p_save(1),p(1)],'r','linewidth',2)

        subplot(2,2,2)
        hold on
        plot([i-1,i],[p_save(2),p(2)],'b','linewidth',2)
        
        subplot(2,2,3)
        hold on
        plot([i-1,i],[p_save(3),p(3)],'b','linewidth',2)

        subplot(2,2,4)
        hold on
        plot([i-1,i],[p_save(4),p(4)],'b','linewidth',2)

        drawnow
        pause(0.2)
    end
    disp('iterazione')
end

% display
disp([p_init; p_true; p])
disp(100*(p-p_true)./p_true)

figure; plot(time,y_pred); hold on; plot(time,glucose,'o')