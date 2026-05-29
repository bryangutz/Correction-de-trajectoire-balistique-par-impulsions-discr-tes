function A = compute_A_numeric(x0, f_dyn)
    n = length(x0);
    A = zeros(n,n);
    eps = 1e-6;

    for i = 1:n
        dx = zeros(n,1);
        dx(i) = eps;

        f_plus  = f_dyn(x0 + dx); % column vector
        f_minus = f_dyn(x0 - dx); % column vector

        A(:,i) = (f_plus - f_minus) / (2*eps);
    end
end