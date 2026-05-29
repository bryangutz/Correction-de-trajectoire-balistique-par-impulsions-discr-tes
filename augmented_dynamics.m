function dz = augmented_dynamics(~, z, n, f_dyn)

    % ---- unpack ----
    x = z(1:n);
    F = reshape(z(n+1:end), n, n);

    % ---- state dynamics ----
    xdot = f_dyn(x);

    % ---- Jacobian ----
    A = compute_A_numeric(x,f_dyn);

    % ---- STM ----
    Fdot = A * F;

    % ---- pack ----
    dz = [xdot; Fdot(:)];
end
