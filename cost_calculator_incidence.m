function cost = cost_calculator_incidence(p, R, p_delta, R_deltas, A_inc)
    [m,n] = size(A_inc);
    r_vec = rot_mat_to_vec(R);
    delta_p_mat = delta_p_cell_to_p_mat(p_delta);
    % create D matrix with size 2mx2n
    D = zeros(2*m, 2*n);
    for k = 1:m
        i = find(A_inc(k,:) == -1);
        D(2*k-1:2*k, 2*i-1:2*i) = -delta_p_mat{k};
    end
    % create p vector from all p's
    p_vec = p(:);
    r_vec = cell2mat(r_vec);
    r_vec = r_vec(:);
    % create A_bar with size 2mx2n
    A_bar = kron(A_inc, eye(2));
    
    U = zeros(2*m, 2*n);
    for k = 1:m
        i = find(A_inc(k,:) == -1);
        j = find(A_inc(k,:) == 1);
        U(2*k-1:2*k, 2*i-1:2*i) = -R_deltas{k};
        U(2*k-1:2*k, 2*j-1:2*j) = eye(2);
    end
    
    cost = norm(A_bar*p_vec + D*r_vec).^2 + norm(U*r_vec).^2;

end