%% Problem 2
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


%2c
% Inner loop roll control
kp_phi = (delta_max_a/e_max_a)*sign(a_phi2);
omega_n_phi = sqrt(abs(a_phi2)*(delta_max_a/e_max_a));
kd_phi = (2*zeta_phi*omega_n_phi - a_phi1)/a_phi2;

% root locus
s = tf('s');
sys = a_phi2 / (s * (s^2 + (a_phi1 + a_phi2 * kd_phi) * s + a_phi2 * kp_phi));
ki_phi = (-100:0.1:100);
figure(4),rlocus(-sys,ki_phi);


%Locus plot, er dette riktig? Skjønner ikke helt overgangen til evans
% men vi ender vel opp med å fjerne det indre integral-leddet til slutt
% uans?
sys = tf([-a_phi2],[1, a_phi1+a_phi2*kd_phi, a_phi2*kp_phi,0]);
figure(1), rlocus(sys); 



% Outer loop control
W_chi = 10; %bandwidht sep, kan endres på (typisk 5-10)
omega_n_chi = omega_n_phi/W_chi; 
zeta_chi = 0.9; %Design parameter
kp_chi = 2*zeta_chi*omega_n_chi*V_g/g; 
ki_chi = (omega_n_chi)^2 * V_g/g;

%2d simulation

sim_time = 500;
d = 1.5*pi/180;
chi_control = timeseries(ones(sim_time,1));
out = sim('model2d',sim_time);
figure(2),plot(out.chi);hold on 
figure(2), plot(out.input);
legend("\chi Out","\chi inp");
title('Course'); ylabel('Angle [deg]');

figure(3),plot(out.delta_a)
title('Aileron'); ylabel('Angle [deg]');
ylim([-50,50]);



