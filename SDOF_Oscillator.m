% Analytical response of damped nonlinear SDOF system - force vibration with IC.
% Comparison between Matlab solvers and own solver. 

clear; clc; dbstop if error; close('all');

%% Definition of input parameters
% System parameters
P.m = 1;                % mass [kg]
P.k = 1;                % stiffness [N/m]
P.c = 0.2;              % damping coefficient [Ns/m]

% Time parameters
Dt = 0.001;             % time step [s]
T = 50;                 % end time [s]
t = 0:Dt:T;             % time [s]
N = size(t,2);          % size of time vector [-]

% Force parameters
P.f = 1;                % magnitude of force [N]
P.Omega = 2;            % radial frequency [rad/s]
F = P.f*cos(P.Omega*t); % force

% Initial conditions
u0 = 1;                 % initial displacement [m]
v0 = 0;                 % initial velocity [m/s]
a0 = (F(1) - P.c*v0 - P.k*u0) / P.m ;    % intial acceleration [m/s^2)
q0 = [u0 v0]';          % state vector
dq0 = [v0 a0]';         % derivative of state vector

%% Prep analytical solution
% Define
wn = sqrt(P.k/P.m);     % Undamped eigenfrequencies
n = P.c/P.m/2;          % Damping
w1 = sqrt(wn^2-n^2);    % Damped eigenfrequencies

% Particular solution amplitudes
Xc = P.f/P.m*(wn^2 - P.Omega^2)/((wn^2 - P.Omega^2)^2 + 4*n^2*P.Omega^2);
Xs = P.f/P.m*2*n*P.Omega/ ((wn^2 - P.Omega^2)^2 + 4*n^2*P.Omega^2);

% General solution amplitudes
A = -Xc + u0;
B = (-n*Xc + n*u0 - Xs*P.Omega + v0)/w1;

%% Get responses
% Total solution function (derivations done in Maple)
W_ana_func_disp = @(t) exp(-n*t).*( A*cos(w1*t) + B*sin(w1*t)) + Xc * cos(P.Omega*t) + Xs*sin(P.Omega*t);
W_ana_func_velo = @(t) ((-A*n+B*w1)*cos(w1*t)-sin(w1*t)*(A*w1+B*n)).* exp(-n*t)-P.Omega*(sin(P.Omega*t)*Xc-Xs*cos(P.Omega*t));
% W_ana_func_acc  = @(t) ((A*n^2-A*w1^2-2*B*n*w1).*cos(w1*t)+...
%     (2*(A*n*w1+(1/2)*B*n^2-(1/2)*B*w1^2)).*sin(w1*t)).*...
%     exp(-n*t)-P.Omega^2*(Xs*sin(P.Omega*t)+Xc*cos(P.Omega*t));

% Compute responses
W_ana_disp = W_ana_func_disp(t);
W_ana_velo = W_ana_func_velo(t);

q = zeros(2,N);
q(:,1) = q0;

for i = 2:N
    q_n = q(:,i-1);
    dq_n = qSDOF(t(i),q_n,P);
    
    q(:,i) = q_n + Dt*dq_n;
end


%% Matlab's solvers
[T,Y] = ode45(@(t_n,q_n) qSDOF(t_n,q_n,P),t,q0);
Y = Y';

%% Plot analytical solutions
figure('units','normalized','outerposition',[0 0.1 0.45 0.6],'PaperPositionMode','auto');
hold on
plot(t, W_ana_disp, 'b','linewidth',1);
plot(t,q(1,:),'r','linewidth',1);
plot(t,Y(1,:),'k-');
hold on
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Comparison of analytical Vs. numerical results')
legend('Analytical','Numerical','ODE45')

%% Velocity: analytical vs numerical 
figure('units','normalized','outerposition',[0 0.1 0.45 0.6],'PaperPositionMode','auto');
hold on
plot(t, W_ana_velo, 'b','linewidth',1);
plot(t,q(2,:),'r','linewidth',1);
plot(t,Y(2,:),'k-');
hold on
xlabel('Time [s]')
ylabel('Velocity [m/s]')
title('Comparison of analytical Vs. numerical results')
legend('Analytical','Numerical','ODE45')

%% Proof of order of FE solver
Dt = logspace(-3,0,7);
T = 50; 
mean_error = zeros(3,length(Dt));

for j = 1:length(Dt)
    t = 0:Dt(j):T;
    N = size(t,2);
    
    q = zeros(2,N);
    q(:,1) = q0;

    for i = 2:N
        q_n = q(:,i-1);
        dq_n = qSDOF(t(i),q_n,P);

        q(:,i) = q_n + Dt(j)*dq_n;
    end
    
    [T1,Y1] = ode45(@(t_n,q_n) qSDOF(t_n,q_n,P),t,q0);
    Y1 = Y1';
    
    [T2,Y2] = ode23(@(t_n,q_n) qSDOF(t_n,q_n,P),t,q0);
    Y2 = Y2';
    
    W_ana_disp = W_ana_func_disp(t);
    W_ana_velo = W_ana_func_velo(t);
    
    error1 = W_ana_disp - q(1,:);
    error2 = W_ana_disp - Y1(1,:);
    error3 = W_ana_disp - Y2(1,:);
    
    mean_error(1,j) = rms(error1);
    mean_error(2,j) = rms(error2);
    mean_error(3,j) = rms(error3);
end

%% Plot & save root mean square error of FE solver vs analytical sol. for different time steps
figure
loglog(Dt,mean_error(1,:),'-o')
xlabel('Time step')
ylabel('Root mean square error')

%% Plot & save RMS error of all solvers vs analytical sol. for different time steps
figure
loglog(Dt,mean_error,'-o')
xlabel('Time step [s]')
ylabel('Root mean square error [m]')
legend('FE solver','ODE45','ODE23')

%% Differences in calculation time
Dt = 0.001;             % time step [s]
T = 50;                 % end time [s]
t = 0:Dt:T;             % time [s]
N = size(t,2);          % size of time vector [-]
% Calculation time of ODE23
t1 = tic;
[T,Y1] = ode23(@(t_n, q_n) qSDOF(t_n, q_n, P), [0 t(end)], q0);
fprintf(1,'ode23 took %.6f [s]\n',toc(t1));
% Calculation time of ODE45
t1 = tic;
[T,Y1] = ode45(@(t_n, q_n) qSDOF(t_n, q_n, P), [0 t(end)], q0);
fprintf(1,'ode45 took %.6f [s]\n',toc(t1));
% Calculation time of own solver
t1 = tic;
q = zeros(2,N);
q(:,1) = q0;
for i = 2:N
    q_n = q(:,i-1);
    dq_n = qSDOF(t(i),q_n,P);
    
    q(:,i) = q_n + Dt*dq_n;
end
fprintf(1,'FE solver took %.6f [s]\n',toc(t1));

% Calculation time analytical
t1 = tic;
W_ana_disp = W_ana_func_disp(t);
W_ana_velo = W_ana_func_velo(t);
fprintf(1,'Analytical calc. took %.6f [s]\n',toc(t1));