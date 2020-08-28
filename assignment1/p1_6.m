clear variable;
% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions. The MSS m-files must be on your
% Matlab path in order to run the script.
%
% System:                      .
%                              q = T(q)w
%                              .
%                            I w - S(Iw)w = tau
% Control law:                        
%                                     ~    ~
%                            tau = -Kdw -kpe 
% 
% Definitions:             
%                            I = inertia matrix (3x3)
%                            S(w) = skew-symmetric matrix (3x3)
%                            T(q) = transformation matrix (4x3)
%                            tau = control input (3x1)
%                            w = angular velocity vector (3x1)
%                            q = unit quaternion vector (4x1)
%
% Author:                   2018-08-15 Thor I. Fossen and Håkon H. Helgesen

%% USER INPUTS
h = 0.1;                     % sample time (s)
N  = 3000;                    % number of samples. Should be adjusted

% model parameters
m = 180;
r = 2;
I = m*r^2*eye(3);            % inertia matrix
I_inv = inv(I);

% controller parameters
kp = 20;
kd = 400;

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

phi = -5*deg2rad;            % initial Euler angles
theta = 10*deg2rad;
psi = -20*deg2rad;

q = euler2q(phi,theta,psi);   % transform initial Euler angles to q

w = [0 0 0]';                 % initial angular rates

table = zeros(N+1,20);        % memory allocation

%% FOR-END LOOP
for i = 1:N+1
   t = (i-1)*h;                  % time
   
   % Control and reference signals
   phi_d = 0*deg2rad;
   theta_d = 15*cos(0.1*t)*deg2rad;
   psi_d = 10*sin(0.05*t)*deg2rad;
   
   dphi_d = 0*deg2rad;
   dtheta_d = -15*0.1*sin(0.1*t)*deg2rad;
   dpsi_d = 10*0.05*cos(0.05*t)*deg2rad;
   
   dTheta_d = [dphi_d dtheta_d dpsi_d]';
   
   trans = Tzyx(phi_d, theta_d);
   w_d = trans\dTheta_d;
   
   w_tilde = w-w_d;
   
   qd = euler2q(phi_d,theta_d,psi_d);
   qd_bar = [qd(1) -qd(2) -qd(3) -qd(4)]';
   
   eta1 = qd_bar(1);
   eta2 = q(1);
   e1 = [qd_bar(2) qd_bar(3) qd_bar(4)]';
   e2 = [q(2) q(3) q(4)]';
   
   eta_tilde = eta1*eta2-e1'*e2;
   e_tilde = eta1*e2+eta2*e1+Smtrx(e1)*e2;
   
   q_tilde = [eta_tilde e_tilde(1) e_tilde(2) e_tilde(3)]';
   
   tau = [-kd*w_tilde(1)-kp*q_tilde(2) -kd*w_tilde(2)-kp*q_tilde(3) -kd*w_tilde(3)-kp*q_tilde(4)]';  % control law

   % Simulation
   [phi,theta,psi] = q2euler(q); % transform q to Euler angles
   [J,J1,J2] = quatern(q);       % kinematic transformation matrices
   
   q_dot = J2*w;                        % quaternion kinematics
   w_dot = I_inv*(Smtrx(I*w)*w + tau);  % rigid-body kinetics
   
   table(i,:) = [t q' phi theta psi w' tau', phi_d, theta_d, psi_d, w_d'];  % store data in table
   
   q = q + h*q_dot;	             % Euler integration
   w = w + h*w_dot;
   
   q  = q/norm(q);               % unit quaternion normalization
end 

%% PLOT FIGURES
t       = table(:,1);  
q       = table(:,2:5); 
phi     = rad2deg*table(:,6);
theta   = rad2deg*table(:,7);
psi     = rad2deg*table(:,8);
w       = rad2deg*table(:,9:11);  
tau     = table(:,12:14);
phi_d   = rad2deg*table(:,15);
theta_d = rad2deg*table(:,16);
psi_d   = rad2deg*table(:,17);
w_d     = rad2deg*table(:,18:20);


figure (1); clf;
hold on;
plot(t, phi, 'b');
plot(t, theta, 'r');
plot(t, psi, 'g');
plot(t, phi_d,'--b');
plot(t, theta_d,'--r');
plot(t, psi_d,'--g');
hold off;
grid on;
legend('\phi', '\theta', '\psi', '\phi_d', '\theta_d', '\psi_d');
title('Euler angles');
xlabel('time [s]'); 
ylabel('angle [deg]');

figure (2); clf;
hold on;
plot(t, w(:,1), 'b');
plot(t, w(:,2), 'r');
plot(t, w(:,3), 'g');
plot(t, w_d(:,1), '--b');
plot(t, w_d(:,2), '--r');
plot(t, w_d(:,3), '--g');
hold off;
grid on;
legend('x', 'y', 'z', 'x_d', 'y_d', 'z_d');
title('Angular velocities');
xlabel('time [s]'); 
ylabel('angular rate [deg/s]');

figure (3); clf;
hold on;
plot(t, tau(:,1), 'b');
plot(t, tau(:,2), 'r');
plot(t, tau(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Control input');
xlabel('time [s]'); 
ylabel('input [Nm]');