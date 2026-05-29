clear; close all; clc;

%% --- Environment constants (from thesis Table 2.2.1) ---
rho0 = 1.225;      % kg/m^3
a0   = 340.429;    % m/s (speed of sound at ground)
Rearth = 6.356766e6; % m
Lat = deg2rad(45);   % latitude used in thesis
g0   = 9.80665*(1-0.0026*cos(2*Lat));    % m/s^2

%% --- Projectile parameters (155 mm from thesis Table 2.3.2) ---
D  = 0.155;            % caliber (m)
S  = 1.89e-2;          % reference area (m^2)
M  = 43.25;            % mass (kg)
Il = 0.15;             % longitudinal moment of inertia (kg m^2)
It = 1.61;             % transverse moment of inertia (kg m^2)

%% --- Example aerodynamic coefficient tables (Mach -> coefficients)
% NOTE: Replace these arrays with the real lookup tables from the thesis.
Mach_table = [0.0 0.8 1.0 1.2 2.0 3.0];
CD_table   = [0.30 0.28 0.46 0.38 0.30 0.28];   % placeholder
CL_table   = [0.5 1.8 1.5 1.7 2.0 1.8]; % placeholder
Cmagf_table= [0.00 -0.08 -0.12 -0.15 -0.13 -0.10]; % placeholder
Cmagm_table= [0.00  0.04  0.06  0.07  0.06  0.05]; % placeholder
CMq_table  = [-22 -20 -18 -15 -12 -10]; % placeholder
CMa_table  = [1.8 2.2 2.6 2.8 2.5 2.2]; % placeholder
Cspin_table= [-0.02 -0.03 -0.04 -0.05 -0.045 -0.040]; % placeholder

% interpolation helper
interpCoeff = @(tableMach, tableVals, Mach) interp1(tableMach, tableVals, Mach, 'linear', 'extrap');

%% --- Initial conditions (example) ---
% Position (local frame L). Use NED-like: x-forward, y-right, z-down (here z is altitude negative)
x0 = 0; y0 = 0; h0 = 0;   % start at ground level (h = altitude)
z0 = -h0;                 % local z coordinate (thesis uses z negative down convention)
V0 = 800;                 % initial speed m/s (gun launch example)
elev = deg2rad(45);       % launch elevation
azim = deg2rad(0);        % azimuth
vx0 = V0 * cos(elev) * cos(azim);
vy0 = V0 * cos(elev) * sin(azim);
vz0 = -V0 * sin(elev);    % negative downwards if z down convention; we'll use consistent sign

% Attitude: use quaternion initial corresponding to yaw, pitch, roll
yaw0 = azim; pitch0 = elev; roll0 = 0;
q0 = eul2quat([yaw0, pitch0, roll0],'ZYX'); % MATLAB builtin: [yaw, pitch, roll]

% Angular rates in body (p,q,r) initial (spin + small transverse)
p0 = 1000; % spin rate rad/s (example — high spin)
q0_ang = 0; r0_ang = 0;

% Full state vector for ODE:
% state = [x; y; z; vx; vy; vz; q0; q1; q2; q3; p; q; r]
state0 = [x0; y0; z0; vx0; vy0; vz0; q0(:); p0; q0_ang; r0_ang];

%% --- Time span ---
T = 75;
tspan = [0 T]; % seconds

%% --- Run simulation ---
[t, X] = ode45(@(t,x) proj_dynamics(t,x, D, S, M, Il, It, ...
                                   rho0, g0, a0, Rearth, Mach_table, CD_table, CL_table, ...
                                   Cmagf_table, Cmagm_table, CMq_table, CMa_table, ...
                                   Cspin_table, interpCoeff), tspan, state0);

%% --- Extract results ---
x = X(:,1); y = X(:,2); z = X(:,3);
vx = X(:,4); vy = X(:,5); vz = X(:,6);
quat = X(:,7:10);
p = X(:,11); q = X(:,12); r = X(:,13);

altitude = -z; % if z is negative down

%% --- Plots ---
figure;
plot3(x, y, altitude); grid on; axis equal;
xlabel('x [m]'); ylabel('y [m]'); zlabel('altitude [m]');
title('3D trajectory');

% figure;
% subplot(3,1,1); plot(t, sqrt(vx.^2+vy.^2+vz.^2)); ylabel('speed [m/s]'); grid on;
% subplot(3,1,2); plot(t, altitude); ylabel('altitude [m]'); grid on;
% subplot(3,1,3); plot(t, rad2deg([p q r])); ylabel('p,q,r [deg/s]'); xlabel('t [s]'); grid on;

figure;
plot(t,sqrt(quat(:,1).^2+quat(:,2).^2+quat(:,3).^2+quat(:,4).^2));ylabel("Norm quaternion"); xlabel("t [s]");
grid on;
%% --- Jacobian A(t) evaluation ---
N = 100; % can be increased

t_table = linspace(0,T,N+1).';
sol = ode45(@(t,x) proj_dynamics(t,x, D, S, M, Il, It, ...
                                   rho0, g0, a0, Rearth, Mach_table, CD_table, CL_table, ...
                                   Cmagf_table, Cmagm_table, CMq_table, CMa_table, ...
                                   Cspin_table, interpCoeff), tspan, state0);
X_table = deval(sol, t_table).'; % row vector

n = size(X_table,2); % # states
A_table = zeros(n,n,N+1);

for k = 1:N+1
    xk = X_table(k,:).';
    A_table(:,:,k) = compute_A_numeric(xk);
end

%% --- Stability analysis ---
Nt = length(t_table);

A_norm = zeros(Nt,1);
real_part = zeros(n,Nt);

for k = 1:Nt
    Ak = A_of_t(t_table(k),t_table,A_table);
    A_norm(k) = norm(Ak,2);
    lambda = eig(Ak);
    real_part(:,k) = real(lambda);
end

figure;
plot(t_table,A_norm,'LineWidth', 1.5);
xlabel('t [s]'); ylabel('||A(t)||_2');
grid on;

figure;
plot(t_table,max(real_part,[],1), 'LineWidth', 1.5);
xlabel('t [s]'); ylabel('max Re(\lambda)');
grid on;

%% --- Forward sensitivity ---
A_handle = @(t) A_of_t(t, t_table, A_table);
F0 = eye(n);
F0vec = F0(:);

[tF,Fvec] = ode45(@(t,F) forward_sensitivity(t,F,A_handle,n),tspan,F0vec);

F_t = zeros(n,n,length(tF));

for k = 1:length(tF)
    F_t(:,:,k) = reshape(Fvec(k,:),n,n);
end

Ft_norm = zeros(length(tF),1);

for k=1:length(tF)
    Ft_norm(k) = norm(F_t(:,:,k),2);
end

figure;
plot(tF,Ft_norm,'LineWidth', 1.5);
xlabel('t [s]'); ylabel('||F(t)||_2');
grid on;
