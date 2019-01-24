function [A, p, R_cells, x_lamb] = incidence_function(opt, i_max)
    % N = 10;ü
    N = opt.N;
    M = opt.M;

    A = zeros(N,N);

    for i = 1:N
        for j = 1:N
            prob = rand(1);
            if(prob < 0.1) && (i ~= j) && (A(i,j) == 0)
    %             A(i,j) = 1;
    %             A(j,i) = -1;
            end
        end
        if(i < N)
            A(i, i+1) = 1;
    %         A(i+1, i) = -1;
        end
    end

    pmax = 10;
    p = pmax*rand(2,N);
    p(:, 1) = [0, 0];
    theta = 2*pi*randn(1,N) - pi;
    R_cells = euler_to_rot_mat(theta);
    R_cells{1} = eye(2);
    [rows, cols] = find(A == 1); M = length(rows);

    A_inc = zeros(M,N);
    for k = 1:M
        i = rows(k); j = cols(k);
        A_inc(k, j) = 1;
        A_inc(k, i) = -1;
    end
    disp(A_inc);

    delta_p_cell = cell(1,M);
    R_delta_cell = cell(1,M);

    for k = 1:M
        i = rows(k); j = cols(k);
        delta_p_cell{k} = R_cells{i}'*(p(:,j) - p(:,i));
        R_delta_cell{k} = R_cells{i}'*R_cells{j};
    end

    delta_p_cell_noise = cell(1,M);
    R_delta_cell_noise = cell(1,M);

    for k = 1:M
        i = rows(k); j = cols(k);
        delta_p_cell_noise{k} = R_cells{i}'*(p(:,j) - p(:,i)) + 0.1*randn(2,1);
        R_noise = euler_to_rot_mat(0.1*(2*pi*randn(1) - pi));
        R_delta_cell_noise{k} = R_cells{i}'*R_cells{j}*R_noise{1};
    end


    % cost = cost_calculator_incidence(p, R_cells, delta_p_cell, R_delta_cell, A_inc)
    % cost_noise = cost_calculator_incidence(p, R_cells, delta_p_cell_noise, R_delta_cell_noise, A_inc)

    [x_vec, W] = cost_calculator_incidence_complex(p, R_cells, delta_p_cell, R_delta_cell, A_inc);
    % cost = real(cost)
    % cost_noise = cost_calculator_incidence_complex(p, R_cells, delta_p_cell_noise, R_delta_cell_noise, A_inc)
    % cost_noise = real(cost_noise)

    [m,n] = size(A_inc);
    
    W_l = @(lambda) W - diag([zeros(n-1,1); lambda]);
%     cost = @(x, lambda) real( x'*W*x + lambda'*(1 - conj(x(n:end)).*x(n:end)));
%     cost_2 = @(x) real( x'*W*x + sum(1 - conj(x(n:end)).*x(n:end)).^2);
%     cost(x_vec, lambda_);

    lambda_0 = ones(n, 1);
    f = @(x) real(-sum(x) );
    f_grad = @(x) gradient_calculater(f, x);
    if(strcmp(opt.method, 'penalty'))
        h = @(x) 0;
    elseif(strcmp(opt.method, 'barrier'))
        h = @(x) inf;
    end
    h_grad = @(x) 0;
    g = @(x) sum(min(eig(W_l(x)), 0).^2);
    g_grad = @(x) gradient_calculater(g,x);
    
    constraints = {h, h_grad, g, g_grad, opt.solver, opt.search};
    
    if(strcmp(opt.method, 'penalty'))
        [lamb_sol, i] = penalty(f, f_grad, constraints, M, lambda_0, i_max);
    elseif(strcmp(opt.method, 'barrier'))
        [lamb_sol, i] = barrier(f, f_grad, constraints, M, lambda_0, i_max);
    end
    
    W_lamb = W_l(lamb_sol);
    [U,S,V] = svd(W_lamb);
    x_lamb = V(:,end);
    if(length(x_lamb) ~= 0)
        for i = n:2*n-1
            x_lamb(i) = x_lamb(i)/norm(x_lamb(i));
        end
    end
    % display(x_lamb);

    % find lambda with cvx sdp
    % cvx_begin
    %     variable lambda(n)
    %     maximize( sum(lambda) )
    %     W_l(lambda) == semidefinite( 2*n-1 )
    %     norm(lambda) <= 100
    % cvx_end

    
end