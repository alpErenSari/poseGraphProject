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
    
    % stack p and r together
    x_vec = [p_vec; r_vec];
    % reduce the x_vec 
    
    % construct the elements of W
    L = A_inc'*A_inc;
    L_bar = kron(L, eye(2));
    Q = D'*D + U'*U;
    W = [L_bar A_bar'*D; D'*A_bar Q];
    
    % construct the W for anchored A
    x_vec_anc = x_vec(3:end);
    A_inc_anc = A_inc(:, 2:end);
    A_inc_anc_bar = kron(A_inc_anc, eye(2));
    L_anc = A_inc_anc_bar'*A_inc_anc_bar;
    S_anc = A_inc_anc_bar'*D;
    W_anc = [L_anc S_anc; S_anc' Q]; 
    
    
%     cost = norm(A_bar*p_vec + D*r_vec).^2 + norm(U*r_vec).^2;
%     cost = x_vec'*W*x_vec;
    cost = x_vec_anc'*W_anc*x_vec_anc;

end