% Project in TTK4190 Guidance and Control of Vehicles 
%
% Author:           Herman Kolstad Jakobsen
%                   Aksel Heggernes
%                   Sondre Holm Fyhn
%                   Iver Myklebust
%
% Study program:    MTTK

clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h  = 0.1;           % sampling time [s]
Ns = 80000;         % no. of samples

sim_with_est = 1;   % decide whether to sim with estimated values or direct noisy measurement (1 or 0)

psi_ref = 0;            % desired yaw angle (rad)
U_d = 7;                % desired cruise speed (m/s)
               
% ship parameters 
m = 17.0677e6;          % mass (kg)
Iz = 2.1732e10;         % yaw moment of inertia about CO (kg m^3)
xg = -3.7;              % CG x-ccordinate (m)
L = 161;                % length (m)
B = 21.8;               % beam (m)
T = 8.9;                % draft (m)
[KT, KQ] = wageningen(0,1.5,0.65,4);               % propeller coefficient (Wageningen) (-)
Dia = 3.3;              % propeller diameter (m)
rho = 1025;             % density of water (kg/m^3)
visc = 1e-6;            % kinematic viscousity at 20 degrees (m/s^2)
eps = 0.001;            % a small number added to ensure that the denominator of Cf is well defined at u=0
k = 0.1;                % form factor giving a viscous correction
t_thr = 0.05;           % thrust deduction number

% rudder limitations
delta_max  = 40 * pi/180;        % max rudder angle      (rad)
Ddelta_max = 5  * pi/180;        % max rudder derivative (rad/s)

% added mass matrix about CO
Xudot = -8.9830e5;
Yvdot = -5.1996e6;
Yrdot =  9.3677e5;
Nvdot =  Yrdot;
Nrdot = -2.4283e10;
MA = -[ Xudot 0    0 
        0 Yvdot Yrdot
        0 Nvdot Nrdot ];

% rigid-body mass matrix
MRB = [ m 0    0 
        0 m    m*xg
        0 m*xg Iz ];
    
Minv = inv(MRB + MA); % Added mass is included to give the total inertia

% ocean current in NED
Vc = 1;                             % current speed (m/s)
betaVc = deg2rad(45);               % current direction (rad)

% wind expressed in NED
Vw = 10;                   % wind speed (m/s)
betaVw = deg2rad(135);     % wind direction (rad)
rho_a = 1.247;             % air density at 10 deg celsius
cy = 0.95;                 % wind coefficient in sway
cn = 0.15;                 % wind coefficient in yaw
A_Lw = 10 * L;             % projected lateral area

% linear damping matrix (only valid for zero speed)
T1 = 20; T2 = 20; T6 = 10;

Xu = -(m - Xudot) / T1;
Yv = -(m - Yvdot) / T2;
Nr = -(Iz - Nrdot)/ T6;
D = diag([-Xu -Yv -Nr]);         % zero speed linear damping

% rudder coefficients (Section 9.5)
b = 2;
AR = 8;
CB = 0.8;

lambda = b^2 / AR;
tR = 0.45 - 0.28*CB;
CN = 6.13*lambda / (lambda + 2.25);
aH = 0.75;
xH = -0.4 * L;
xR = -0.5 * L;

X_delta2 = 0.5 * (1 - tR) * rho * AR * CN;
Y_delta = 0.25 * (1 + aH) * rho * AR * CN; 
N_delta = 0.25 * (xR + aH*xH) * rho * AR * CN;   

% input matrix
Bu = @(u_r,delta) [ (1-t_thr)  -u_r^2 * X_delta2 * delta
                        0      -u_r^2 * Y_delta
                        0      -u_r^2 * N_delta            ];
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
% Heading Controller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rudder control law
if sim_with_est == 1
   wb = 0.03;
   zeta = 0.2;  
else
   wb = 0.06;
   zeta = 1; 
end
wn = 1 / sqrt( 1 - 2*zeta^2 + sqrt( 4*zeta^4 - 4*zeta^2 + 2) ) * wb;

% linearized sway-yaw model (see (7.15)-(7.19) in Fossen (2021)) used
% for controller design. The code below should be modified.
M_lin = MRB + MA;
M_lin = M_lin(2:3,2:3);
M_lin_inv = inv(M_lin);

CRB_lin = [ 0    m*U_d
            0    m*xg*U_d  ];
        
CA_lin = [  0                       -Xudot*U_d
            -Yvdot*U_d+Xudot*U_d    -Yrdot*U_d  ];
        
N_lin = CRB_lin + CA_lin + D(2:3,2:3);
b_lin = [-2*U_d*Y_delta -2*U_d*N_delta]';

[num, den] = ss2tf(-M_lin_inv*N_lin, M_lin_inv*b_lin, [0 1], 0);

poles_lin = roots(den);
zeros_lin = roots(num);

T3_lin = -1/zeros_lin;   
T1_lin = -1/poles_lin(1);
T2_lin = -1/poles_lin(2);

T_lin = T1_lin+T2_lin-T3_lin;       % nomoto first-order time constant, eq. (7.24)
K_lin = num(3)/den(3);              % gain can be found by using the steady-state value (s=0) of the transfer function

% controller gains (example 15.7 / Algorithm 15.1)
Kp = wn^2*T_lin/K_lin;         
Kd = (2*zeta*wn*T_lin-1)/K_lin;
Ki = wn/10*Kp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
% State Estimation (Kalman Filter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system matrices
A_est = [   0   1           0 
            0   -1/T_lin    -K_lin/T_lin
            0   0           0            ];
        
B_est = [0  K_lin/T_lin     0]';

E_est = [   0   0           
            1   0    
            0   1   ];
        
C_est = [1  0   0];

%disp(rank(obsv(A_est, C_est)));

% discretized matrices (first order)
A_d = eye(3) + h*A_est;
B_d = h*B_est;
C_d = C_est;
D_d = 0;
E_d = h*E_est;

% standard deviation for measurement noise
sigma_psi = deg2rad(0.5);
sigma_r = deg2rad(0.1);

% covariance matrices for process and measurement noise
sigma_q1 = 1 / 10000;   % should be tuned to get desired convergence
sigma_q2 = 1 / 100000;   % should be tuned to get desired convergence

Q_d = diag([sigma_q1^2 sigma_q2^2]);
R_d = sigma_psi^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
% Initial states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta = [0 0 0]';
nu  = [0.1 0 0]';
nu_c = [0 0 0]';    % current velocities
xd = [0 0 0 ]';     % initial reference model states (3rd order)
z = 0;              % integral state (PID controller)
delta = 0;
n = 0;
Qm = 0;             % produced torque by main motor (Nm)
wp = 2;             % which waypoint the ship is currently targeting for
WP = load('WP.mat'); % waypoints
WP = WP.WP;
chi_d = 0; 
y_int = 0; 
delta_los = 1000;
kappa = 4;

% initial states for Kalman filter
x0 = [0 0 0]';                 % [yaw_angle, yaw_rate, rudder_bias]
P0 = diag([1/1000 1/1000000 1/10000]);     % initial covariance matrix

x_prd = x0;
P_prd = P0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(Ns+1,24);                % table of simulation data

for i=1:Ns+1
    % Noisy measurement
    psi_meas = eta(3) + normrnd(0,sigma_psi);
    r_meas = nu(3) + normrnd(0,sigma_r);
    
    eta(3) = wrapTo2Pi(eta(3));         % solve "plotting bug" of desired vs actual yaw angle
    
    t = (i-1) * h;                      % time (s)
    R = Rzyx(0,0,eta(3));
    
    % current (should be added here)
    nu_c(1) = Vc*cos(betaVc - eta(3));  % surge current
    nu_c(2) = Vc*sin(betaVc - eta(3));  % sway current
    
    nu_r = nu - nu_c;
    
    % wind (should be added here)
    if t > 200
        gamma_w = eta(3) - betaVw - pi;
        C_Ywind = cy*sin(gamma_w);              % sway wind coefficient
        C_Nwind = cn*sin(2*gamma_w);            % yaw wind coefficient
        Ywind = 1/2*rho_a*Vw^2*C_Ywind*A_Lw;    % expression for wind moment in sway should be added.
        Nwind = 1/2*rho_a*Vw^2*C_Nwind*A_Lw*L;  % expression for wind moment in yaw should be added.
    else
        Ywind = 0;
        Nwind = 0;
    end
    tau_env = [0 Ywind Nwind]';
    
    % state-dependent time-varying matrices
    CRB = m * nu(3) * [ 0 -1 -xg 
                        1  0  0 
                        xg 0  0  ];
                    
    % coriolis due to added mass
    CA = [  0   0   Yvdot * nu_r(2) + Yrdot * nu_r(3)
            0   0   -Xudot * nu_r(1) 
          -Yvdot * nu_r(2) - Yrdot * nu_r(3)    Xudot * nu_r(1)   0];
    N = CRB + CA + D;
    
    % nonlinear surge damping
    Rn = L/visc * abs(nu_r(1));
    Cf = 0.075 / ( (log(Rn) - 2)^2 + eps);
    Xns = -0.5 * rho * (B*L) * (1 + k) * Cf * abs(nu_r(1)) * nu_r(1);
    
    % cross-flow drag
    Ycf = 0;
    Ncf = 0;
    dx = L/10;
    Cd_2D = Hoerner(B,T);
    for xL = -L/2:dx:L/2
        vr = nu_r(2);
        r = nu_r(3);
        Ucf = abs(vr + xL * r) * (vr + xL * r);
        Ycf = Ycf - 0.5 * rho * T * Cd_2D * Ucf * dx;
        Ncf = Ncf - 0.5 * rho * T * Cd_2D * xL * Ucf * dx;
    end
    d = -[Xns Ycf Ncf]';
    
    % guidance
    [~,number_wp] = size(WP);
    
    if wp <= number_wp

        wp_x1 = WP(1, wp-1);
        wp_y1 = WP(2, wp-1);
        
        wp_x2 = WP(1, wp);
        wp_y2 = WP(2, wp);
        ship_x = eta(1);
        ship_y = eta(2);
        
        % ILOS Guidance
        [chi_d, y_e] = guidanceInt(wp_x1, wp_y1, wp_x2, wp_y2, ship_x, ship_y, delta_los, kappa, y_int); 
        y_int_dot = (delta_los * y_e) / (delta_los^2 + (y_e+kappa*y_int)^2);
        y_int = euler2(y_int_dot,y_int,h);
        
        dist_to_wp = norm([wp_x2, wp_y2] - [ship_x ship_y]);
        if dist_to_wp < 2000
            wp = wp + 1;
        end
    end
    beta_c = asin(nu(2)/norm([nu(1) nu(2)]));   % crab angle
    beta = asin((nu(2)-nu_c(2))/norm([nu(1)-nu_c(1) nu(2)-nu_c(2)])); % sideslip
    
    chi = eta(3) + beta_c;                      % course
    psi_ref = chi_d;                            % Setting psi ref
    
    % 3rd-order reference model for yaw, eq.(12.12)
    wref = 0.05;    % natural frequency for reference model
    Ad = [ 0 1 0
           0 0 1
           -wref^3  -3*wref^2  -3*wref ];
    Bd = [0 0 wref^3 ]';
    xd_dot = Ad * xd + Bd * psi_ref;
    
    psi_d = xd(1);  % desired yaw angle (rad)
    r_d = xd(2);    % desired yaw rate (rad/s)
    u_d = U_d;      % desired cruise spees (m/s)
    
    % thrust 
    thr = rho * Dia^4 * KT * abs(n) * n;    % thrust command (N)
    
    % torque
    Q = rho * Dia^4 * KQ * abs(n) * n;      % torque command (Nm)
    
    % kalman filter
    K = P_prd * C_d' * inv(C_d * P_prd * C_d' + R_d);   % filter gain
    IKC = eye(size(x0,1)) - K* C_d;
    
    y = psi_meas;   % measurement
    
    x_hat = x_prd + K * ( ssa(y - C_d * x_prd) - D_d * delta ); % corrector
    P_hat = IKC * P_prd * IKC' + K * R_d * K';
    
    x_prd = A_d * x_hat + B_d * delta;      % predictor
    P_prd = A_d * P_hat * A_d' + E_d * Q_d * E_d';
    
    x_hat(1) = wrapTo2Pi(x_hat(1));        % fix plotting bug
        
    % control law
    if sim_with_est == 1
       e_psi = ssa(x_hat(1) - psi_d);  % estimated value
       e_r = x_hat(2) - r_d;
    else
       e_psi = ssa(psi_meas-psi_d);    % direct noisy measurement
       e_r = r_meas-r_d;  
    end
    
    delta_c_unsat = -Kp*e_psi-Ki*z-Kd*e_r;    % rudder angle command (rad)
    
    % ship dynamics
    u = [ thr delta ]';
    tau = Bu(nu_r(1),delta) * u;
    nu_dot = [nu(3)*nu_c(2) -nu(3)*nu_c(1) 0]' + Minv * (tau_env + tau - N * nu_r - d); 
    eta_dot = R * nu;    
    
    % Rudder saturation and dynamics (Sections 9.5.2)
    if abs(delta_c_unsat) >= delta_max
        delta_c = sign(delta_c_unsat)*delta_max;
        z = z-(h/Ki)*(delta_c-delta_c_unsat);       % anti-windup
    else
        delta_c = delta_c_unsat;
    end
    
    delta_dot = delta_c - delta;
    if abs(delta_dot) >= Ddelta_max
        delta_dot = sign(delta_dot)*Ddelta_max;
    end    
    
    % propeller dynamics
    Im = 100000; Tm = 10; Km = 0.6;             % propulsion parameters
    
    % added feedforward
    Td = (u_d-nu_c(1))*Xu/(t_thr-1);            % desired thrust (N)
    
    n_term = Td/(rho * Dia^4 * KT);             
    n_d = sign(n_term) * sqrt(abs(n_term));     % desired propeller speed (rps)
                                                
    Qf = 0;                                     % friction torque (Nm)
    Qd =  rho * Dia^4 * KQ * abs(n_d) * n_d;    % desired propeller moment (Nm)
    Y = Qd/Km;                                  % control input to main motor
    
    Qm_dot = -Qm/Tm + Km/Tm*Y;
    n_dot = (Qm-Q-Qf)/Im;                      
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [t n_d delta_c n delta eta' nu' u_d psi_d r_d z beta_c beta chi chi_d psi_meas r_meas x_hat'];       
     
    % Euler integration
    xd = euler2(xd_dot,xd,h);               % reference model
    z = euler2(e_psi,z,h);                  % integral state
    Qm = euler2(Qm_dot,Qm,h);
    eta = euler2(eta_dot,eta,h);
    nu  = euler2(nu_dot,nu,h);
    delta = euler2(delta_dot,delta,h);   
    n  = euler2(n_dot,n,h);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t       = simdata(:,1);                 % s
n_d     = 60 * simdata(:,2);            % rpm
delta_c = (180/pi) * simdata(:,3);      % deg
n       = 60 * simdata(:,4);            % rpm
delta   = (180/pi) * simdata(:,5);      % deg
x       = simdata(:,6);                 % m
y       = simdata(:,7);                 % m
psi     = (180/pi) * simdata(:,8);      % deg
u       = simdata(:,9);                 % m/s
v       = simdata(:,10);                % m/s
r       = (180/pi) * simdata(:,11);     % deg/s
u_d     = simdata(:,12);                % m/s
psi_d   = (180/pi) * simdata(:,13);     % deg
r_d     = (180/pi) * simdata(:,14);     % deg/s
z       = simdata(:,15); 
beta    = (180/pi) * simdata(:,16);     % deg
beta_c  = (180/pi) * simdata(:,17);     % deg
chi     = (180/pi) * simdata(:,18);     % deg
chi_d   = (180/pi) * simdata(:,19);     % deg
psi_meas = (180/pi) * simdata(:,20);    % deg
r_meas   = (180/pi) * simdata(:,21);    % deg/s
psi_est = (180/pi) * simdata(:,22);     % deg
r_est   = (180/pi) * simdata(:,23);     % deg/s
bias_est = (180/pi) * simdata(:,24);    % deg


figure(1)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions (m)'); xlabel('time (s)'); 
subplot(312)
plot(t,psi,t,psi_d,'linewidth',2);
title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
subplot(313)
plot(t,r,t,r_d,'linewidth',2);
title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocities (m/s)'); xlabel('time (s)');
subplot(312)
plot(t,n,t,n_d,'linewidth',2);
title('Actual and commanded propeller speed (rpm)'); xlabel('time (s)');
subplot(313)
plot(t,delta,t,delta_c,'linewidth',2);
title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');

figure(3) 
figure(gcf)
subplot(211)
plot(t,u,'linewidth',2);
title('Actual surge velocity (m/s)'); xlabel('time (s)');
subplot(212)
plot(t,v,'linewidth',2);
title('Actual sway velocity (m/s)'); xlabel('time (s)');

figure(4)
figure(gcf)
hold on
siz=size(WP);
for ii=1:(siz(2)-1)   
plot([WP(2,ii), WP(2,ii+1)], [WP(1,ii), WP(1,ii+1)], 'r-x')
end
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions (m)');

% noisy measurement vs true value
% figure(5)
% hold on;
% xlabel('time (s)');
% ylabel('(deg)');
% plot(t, psi_meas);
% plot(t, psi);
% legend('measured', 'true');
% title('True vs measured yaw angle');
% 
% figure(6)
% hold on;
% xlabel('time (s)');
% ylabel('(deg/s)');
% plot(t, r_meas);
% plot(t, r);
% legend('measured', 'true');

% estimated states vs true value 
% figure(7)
% hold on;
% xlabel('time (s)');
% ylabel('(deg)');
% plot(t, psi-psi_est);
% title('Yaw angle error between estimated and true value');
% 
% figure(8)
% hold on;
% xlabel('time (s)');
% ylabel('(deg/s)');
% plot(t, r-r_est);
% title('Yaw rate error between estimated and true value');
% 
% figure(9)
% hold on;
% xlabel('time (s)');
% ylabel('(deg)');
% plot(t, bias_est);
% title('Estimated rudder bias');