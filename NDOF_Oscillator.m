clear; clc; dbstop if error; close('all');

%% Load data
load('NDOF_data.mat')
N = P.N;
d = sqrt(eig(P.K,P.M));
% Display smallest eigenfrequencies
disp(min(d))
disp(d(end-1))

%% Solve problem
[~,W_beam] = ode23(@(t_n,q_n) qNDOF(t_n,q_n,P), t_vec, zeros(2*(N-1),1));
W_beam = [fliplr(W_beam(:,1:N-1)), fliplr(W_beam(:,N:end))];

%% Plot results
% Open figure
figure('units','normalized','outerposition',[0 0.1 0.45 0.85],'PaperPositionMode','auto');

% Plot
title('Cantilever beam subjected to dynamic load')
x = linspace(0,P.L,P.N);
surface(t_vec, x(2:end), rot90(W_beam(:,1:P.N-1)),'edgecolor','none')
xlabel('time [s]')
ylabel('x-loc [m]')
zlabel('displacement [m]')
legend
