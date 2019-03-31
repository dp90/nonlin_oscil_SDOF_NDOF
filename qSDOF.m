function q_n_dot = qSDOF( t_n, q_n, P )

u = q_n(1);
du = q_n(2);
ddu = (P.f*cos(P.Omega*t_n) - P.c*du - P.k*u) / P.m;
q_n_dot = [du ddu]';

end