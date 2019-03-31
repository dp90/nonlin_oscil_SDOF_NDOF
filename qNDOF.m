function q_n_dot = qNDOF( t_n, q_n, P )

% Create vector u, containing displacements at time n
u = q_n(1:length(q_n)/2);

% Create vector v, containing velocities at time n
v = q_n(length(q_n)/2+1:end);

Force = ones(20,1)*cos(2.5*pi*t_n);

% Calculate vector a, containing accelerations at time n (from u, v, M, C, K, and f) 
a = mldivide(P.M,(Force - P.C*v - P.K*u));

% Combine velocity and acceleration vectors in derivative of state vector
q_n_dot = [v; a];

end