clear variables;

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

delta_max = 30*deg2rad;
e_phimax = 15*deg2rad;
zeta_phi = 0.707;

a_phi1 = 2.87;
a_phi2 = -0.65;

wn_phi = sqrt(abs(a_phi2) * delta_max / e_phimax);

% gains
kp_phi = delta_max / e_phimax * sign(a_phi2);
kd_phi = (2 * zeta_phi - a_phi1) / a_phi2;

% root locus
s = tf('s');
sys = a_phi2 / (s * (s^2 + (a_phi1 + a_phi2 * kd_phi) * s + a_phi2 * kp_phi));
ki_phi = (-100:0.1:100);
rlocus(-sys,ki_phi);