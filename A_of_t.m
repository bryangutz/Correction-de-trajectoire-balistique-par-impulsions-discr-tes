function A = A_of_t(tq, t_table, A_table)
    n = size(A_table,1);
    A = zeros(n,n);

    for i = 1:n
        for j = 1:n
            A(i,j) = interp1(t_table, squeeze(A_table(i,j,:)), ...
                             tq, 'linear', 'extrap');
        end
    end
end