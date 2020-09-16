%% Problem 2f
clc;close all; clear;

%Parameters
delta_max_a = 30;
e_max_a = 15; 
zeta_phi = 0.707; 
V_a = 580;
V_g = V_a;
g = 9.81;


% From 2a
a_phi1 = 2.87;
a_phi2 = -0.65;

%From 2c
%inner loop
kp_phi = (delta_max_a/e_max_a)*sign(a_phi2);
omega_n_phi = sqrt(abs(a_phi2)*(delta_max_a/e_max_a));
kd_phi = (2*zeta_phi*omega_n_phi - a_phi1)/a_phi2;

% Outer loop control
W_chi = 10; %bandwidht sep, kan endres på (typisk 5-10)
omega_n_chi = omega_n_phi/W_chi; 
zeta_chi = 0.9; %Design parameter
kp_chi = 2*zeta_chi*omega_n_chi*V_g/g; 
ki_chi = (omega_n_chi)^2 * V_g/g;


%2f simulation

%System matrices
A = [-0.322, 0.052, 0.028, -1.12, 0.002;
             0, 0, 1, -0.001, 0;
             -10.6, 0, -2.87, 0.46, -0.65;
             6.87, 0, -0.04, -0.32, -0.02;
             0, 0, 0, 0, -7.5];
 
B = [0; 0; 0; 0; 7.5;];
C = [ 1 0 0 0 0
      0 1 0 0 0
      0 0 1 0 0
      0 0 0 1 0 ];
  
% Actuator dynamics
H_l = tf([0 7.5],[1 7.5]);

sim_time = 500;
chi_control = timeseries(ones(sim_time,1));
out = sim('model2f',sim_time);
figure(2),plot(out.control);hold on 
figure(2), plot(out.input);
legend("\chi Out","\chi inp");
title('Course'); ylabel('Angle [deg]');

figure(3),plot(out.delta_a)
title('Aileron'); ylabel('Angle [deg]');
ylim([-50,50]);