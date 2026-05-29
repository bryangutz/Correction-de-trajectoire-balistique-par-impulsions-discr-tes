function dFvec = forward_sensitivity(t, Fvec, A_handle, n)
    % reconstruct F(t) as matrix
    F = reshape(Fvec, n, n);

    % Jacobian
    A = A_handle(t);

    Fdot = A * F; 
    dFvec = Fdot(:);
end
