function [x_vec_anc, W_anc] = cost_calculator_incidence_complex(p, R, p_delta, R_deltas, A_inc)
    [m,n] = size(A_inc);
    r_vec = rot_mat_to_vec(R);
    delta_p_mat = delta_p_cell_to_p_mat(p_delta);
    % create D matrix with size 2mx2n
    D = zeros(m, n);
    for k = 1:m
        i = find(A_inc(k,:) == -1);
%         theta_ = atan(delta_p_mat{k}(2,1) / delta_p_mat{k}(1,1));
%         alpha_ = norm(delta_p_mat{k}(:,1));
        D(k, i) = -complex(delta_p_mat{k}(1,1), delta_p_mat{k}(2,1));
    end
    % create p vector from all p's
    p_vec = complex(p(1,:), p(2,:));
    r_vec = cell2mat(r_vec);
    r_vec = complex(r_vec(1,:), r_vec(2,:));
    
    U = zeros(m, n);
    for k = 1:m
        i = find(A_inc(k,:) == -1);
        j = find(A_inc(k,:) == 1);
%         theta_ = atan(R_deltas{k}(2,1) / R_deltas{k}(1,1));
        U(k, i) = -complex(R_deltas{k}(1,1), R_deltas{k}(2,1));
        U(k, j) = 1;
    end
    
    % stack p and r together
    x_vec = [transpose(p_vec); transpose(r_vec)];
    % reduce the x_vec 
    
    % construct the elements of W
    L = A_inc'*A_inc;
    Q = D'*D + U'*U;
    
    % construct the W for anchored A
    x_vec_anc = x_vec(2:end);
    A_inc_anc = A_inc(:, 2:end);
    L_anc = A_inc_anc'*A_inc_anc;
    S_anc = A_inc_anc'*D;
    W_anc = [L_anc S_anc; S_anc' Q]; 
    
    
%     cost = norm(A_bar*p_vec + D*r_vec).^2 + norm(U*r_vec).^2;
%     cost = x_vec'*W*x_vec;
    cost = x_vec_anc'*W_anc*x_vec_anc;
%     cost = 0;

end