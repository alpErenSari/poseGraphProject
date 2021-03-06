function [p, R_matrices, delta_p_cell, R_delta_cell] = random_pose_loop_closure(N, pmax)
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
    delta_p_cell = cell(1,N);
    R_delta_cell = cell(1,N);

    for i = 1:N
    %     R_i = R(1:2, 2*i-1:2*i);
        loop_prob = rand(1);
%         loop_prob = 1;
        
        delta_p = pmax*rand(2,1);
        delta_theta = 2*pi*rand(1) - pi;
        R_delta = euler_to_rot_mat(delta_theta);
        R_i = R_matrices{i};
        R_j = R_i*R_delta{1};
        R_matrices{i+1} = R_j;
        if(loop_prob < 0.1)
            c = randi(i-1);
            p{i+1} = p{i-c};
            delta_p = R_i'*(p{i+1} - p{i});
%             disp('First part');
        else
            p{i+1} = p{i} + R_i*delta_p;
        end
        delta_p_cell{i} = delta_p;
        R_delta_cell{i} = R_delta{1};
        
    end
end