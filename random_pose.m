function [p, R_matrices] = random_pose(N, delta_p, R_deltas)
%     % generate p's between 0 and 10
%     delta_p = p_max*rand(2, N);
%     % generate thetas btw -pi and pi
%     delta_theta = 2*pi*rand(1,N) - pi;
%     R_deltas = euler_to_rot_mat(delta_theta);

    p = cell(1,N+1);
    % R = zeros(2,2*(N+1));
    R_matrices = cell(1,N+1);
    % R(:,1:2) = eye(2);
    p{1} = [0 0]';
    R_matrices{1} = eye(2);

    for i = 1:N
    %     R_i = R(1:2, 2*i-1:2*i);
        R_i = R_matrices{i};
    %     keyboard;
        R_ij = R_deltas{i};
        p{i+1} = p{i} + R_i*delta_p(:, i);
        R_j = R_i*R_ij;
    %     R(:, 2*i+1:2*i+2) = R_j;
        R_matrices{i+1} = R_j;
    end