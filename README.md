# Nonlinear oscillators with a single or n degrees of freedom
Solve state-space vector for nonlinear damped oscillators with a single or n d.o.f.s in time domain. 
Includes comparison of different solvers for SDOF system. 

SDOF_Oscillator.m computes the response of an SDOF system and uses function qSDOF to determine the space-state vector.
NDOF_Oscillator.m computes the response of an NDOF system and uses function qNDOF to determine the space-state vector. It also uses NDOF_data for the eigenfrequencies of the NDOF system. 
