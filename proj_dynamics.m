
%% --- Dynamics function ---
function xdot = proj_dynamics(~, x, D, S, M, Il, It, rho0, g0, a0, Rearth, ...
                             Mach_table, CD_table, CL_table, ...
                             Cmagf_table, Cmagm_table, CMq_table, CMa_table, ...
                             Cspin_table, interpCoeff)
    % Unpack state
    pos = x(1:3);    % [x;y;z]
    Vloc = x(4:6);   % local-frame velocity vector
    q = x(7:10);     % quaternion (scalar-first is what eul2quat gives)
    % q = q/norm(q); % renormalization
    p = x(11); q_ang = x(12); r = x(13); % body rates

    % Recompute altitude and local variables
    z = pos(3); altitude = -z;
    v_rel = Vloc; v = norm(v_rel);
    if v < 1e-6, v = 1e-6; end

    % Air properties at current altitude (simple ISA approximation)
    T0 = 288.16; rho = rho0 * ( (T0 - 0.0065*altitude)/T0 )^4.2561;
    vsound = a0 * sqrt((T0 - 0.0065*altitude)/T0);
    Mach = v / max(vsound,1e-6);

    % Extract body axis 1B (unit vector) from quaternion (R body to local)
    Rbl = quat2rotm(q'); % quaternion->rotation matrix (body to local)
    % 1B is first column of Rbl expressed in local frame.
    e1B_local = Rbl(:,1);

    % incidence angles alpha, beta (attack and sideslip)
    % alpha: rotation about 2B between 1B and velocity vector; compute using vector math:
    Vb = Rbl' * v_rel; % express velocity in body frame
    Vb_hat = Vb / (norm(Vb)+1e-9);
    alpha = atan2( Vb(3), Vb(1) ); % approx: angle of attack (small-angle assumption)
    beta  = asin( Vb(2) / max(norm(Vb),1e-9) ); % sideslip approx
    alpha_t = acos(cos(alpha)*cos(beta)); % total angle of attack

    % Interpolate aerodynamic coefficients from Mach (placeholder tables)
    CD = interpCoeff(Mach_table, CD_table, Mach);
    CL = interpCoeff(Mach_table, CL_table, Mach);
    Cmagf = interpCoeff(Mach_table, Cmagf_table, Mach);
    Cmagm = interpCoeff(Mach_table, Cmagm_table, Mach);
    CMq = interpCoeff(Mach_table, CMq_table, Mach);
    CMa = interpCoeff(Mach_table, CMa_table, Mach);
    Cspin = interpCoeff(Mach_table, Cspin_table, Mach);

    % --- Forces (Table 2.4.1 style) ---
    % Drag: -0.5 * rho * S * CD * v^2* Vhat  (vector opposite to velocity)
    % Vhat = v_rel / v;
    Fdrag = -0.5 * rho * S * CD * v * v_rel;

    % Lift (vector): 0.5*rho*S*CL*( v x (1B x v) )
    % Compute 1B in local frame (e1B_local), use v_rel in local
    cross1 = cross(1*e1B_local, v_rel);
    Flift = 0.5 * rho * S * CL * cross(v_rel, cross1);

    % Magnus force: 0.5*rho*S*(p*D/v)*Cmagf * (v x 1B)
    % Here p is spin (body longitudinal), but p is in body rates (rad/s). We use absolute spin rate p.
    % Express spin axis 1B in local frame
    Fmag = 0.5 * rho * S * (p * D / max(v,1e-6)) * Cmagf * cross(v_rel, e1B_local);

    % Gravity (local frame): mass * g (downwards)
    g = g0 * ( Rearth/(Rearth+altitude) )^2; % approx
    Fgrav = [0;0;M*g];

    % Sum forces in local frame
    Ftot_local = Fdrag + Flift + Fmag + Fgrav;

    % Translational acceleration (approx neglecting Coriolis / Earth rotation here; could add)
    acc_local = Ftot_local / M;

    % --- Moments (body frame) using Table 2.4.2 (expressed in body frame) ---
    % We'll compute moments in local frame and then express in body frame.
    % Use simplified expressions from thesis:
    % Magnus moment: 0.5*rho*S*D*(p*D/v)*Cmagm * (1B x (v x 1B))
    Mmag_local = 0.5 * rho * S * D * (p * D / max(v,1e-6)) * Cmagm * v * cross(e1B_local, cross(v_rel, e1B_local));
    % Overturning: 0.5*rho*S*D*CMa * v * (v x 1B)
    Mover_local = 0.5 * rho * S * D * CMa * v * cross(v_rel, e1B_local);
    % Pitch damping moment approx: 0.5*rho*S*D^2*CMq * v * (1B x (Omega x 1B))
    Omega_body = [p; q_ang; r];
    % convert 1B to body coordinates (it's [1;0;0] in body)
    Mpitch_body = 0.5 * rho * S * D^2 * CMq * v * cross([1;0;0], cross(Omega_body, [1;0;0]));
    % roll damping moment (example)
    Mroll_body = 0.5 * rho * S * D * (p * D / max(v,1e-6)) * Cspin * v^2 * [1;0;0];

    % Sum moments in body frame: convert local moments to body frame
    Mmag_body = Rbl' * Mmag_local;
    Mover_body = Rbl' * Mover_local;
    Mtot_body = Mmag_body + Mover_body + Mpitch_body + Mroll_body;

    % --- Rotational dynamics (Euler equations) ---
    % Using: I * domega/dt + omega x (I*omega) = Mtot
    I = diag([Il, It, It]);
    omega = Omega_body;
    domega = I \ (Mtot_body - cross(omega, I*omega));

    % Assign derivatives
    pos_dot = Vloc;
    vel_dot = acc_local; % local frame acceleration

    % Quaternion kinematics: q_dot = 0.5 * Omega_quat * q
    omega_quat = [0; omega];        % pure quaternion (0, p, q, r)

    lambda = 5; % 1-10
    e = 1-(q.'*q); % error
    q_dot = 0.5 * quatmultiply(q', omega_quat')' + lambda*e*q;

    % Pack derivatives
    xdot = zeros(13,1);
    xdot(1:3) = pos_dot;
    xdot(4:6) = vel_dot;
    xdot(7:10) = q_dot;
    xdot(11:13) = domega;
end
