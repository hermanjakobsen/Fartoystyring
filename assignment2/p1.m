%% System

A = [
    -0.322  0.052   0.028   -1.12   0.002;
    0       0       1       -0.001  0;
    -10.6   0       -2.87   0.46    -0.65;
    6.87    0       -0.04   -0.32   -0.02;
    0       0       0       0       -7.5
];

V_a = 580; 


%% Problem 1c
Y_v = A(1,1); 
Y_r = A(1,4)*V_a;
N_r = A(4,4);
N_v = A(4,1)/V_a;

p_dr_1   = (Y_v + N_r)/2 + sqrt(((Y_v + N_r)/2)^2 - (Y_v*N_r - N_v*Y_r));
p_dr_2   = (Y_v + N_r)/2 - sqrt(((Y_v + N_r)/2)^2 - (Y_v*N_r - N_v*Y_r));

w_0 = sqrt(real(p_dr_1)^2 + imag(p_dr_2)^2);
zeta = -real(p_dr_1)/w_0;


%% Problem 1d
L_v = A(3,1)/V_a; 
L_r = A(3,4);

p_s = (N_r*L_v - N_v*L_r)/L_v;


%% Results check
% Compare results with this output?
damp(A)