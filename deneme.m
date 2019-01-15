clear;clc;

N = 10;
sigma_p = 1;
sigma_r = 0.3;
% generate p's between 0 and 10
delta_p = 10*rand(2, N);
% generate thetas btw -pi and pi
delta_theta = 2*pi*rand(1,N) - pi;
R_deltas = euler_to_rot_mat(delta_theta);

delta_p_noised = delta_p + sigma_p*randn(2,N);
delta_theta_noised = delta_theta + sigma_r*randn(1,N);
R_deltas_noised = euler_to_rot_mat(delta_theta_noised);

% 
% p = cell(1,N+1);
% % R = zeros(2,2*(N+1));
% R_matrices = cell(1,N+1);
% % R(:,1:2) = eye(2);
% p{1} = [0 0]';
% R_matrices{1} = eye(2);
% 
% for i = 1:N
% %     R_i = R(1:2, 2*i-1:2*i);
%     R_i = R_matrices{i};
% %     keyboard;
%     R_ij = R_deltas{i};
%     p{i+1} = p{i} + R_i*delta_p(:, i);
%     R_j = R_i*R_ij;
% %     R(:, 2*i+1:2*i+2) = R_j;
%     R_matrices{i+1} = R_j;
% end

% [p, R_matrices] = random_pose(N, delta_p, R_deltas);
[p, R_matrices, delta_p_cell, R_delta_cell] = random_pose_loop_closure(N, 10);
[p_noised, R_matrices_noised] = random_pose(N, delta_p_noised, R_deltas_noised);

delta_p_mat = cell2mat(delta_p_cell);
cost_org = cost_calculator(p, R_matrices, delta_p_mat, R_delta_cell)
% cost_noised = cost_calculator(p_noised, R_matrices_noised, delta_p_noised, R_deltas_noised)

% figure;
% % plot(p(1,:), p(2,:));
% for i = 1:length(p)
%     plot(p{i}(1), p{i}(2), 'x');
%     text(p{i}(1),p{i}(2),num2str(i));
% %     i=i+1;
%     hold on;
% end
    
s = 1:length(p)-1;
t = 2:length(p);
G = digraph(s,t);

G.Nodes.translation = p';
G.Nodes.rotation = R_matrices';

p_mat = cell2mat(p);
p_mat_noised = cell2mat(p_noised);

figure;
plot(G, 'XData', p_mat(1,:), 'YData', p_mat(2,:));
title('Generated Random Pose Graph');
xlabel('x'); ylabel('y');

% figure;
% plot(G, 'XData', p_mat_noised(1,:), 'YData', p_mat_noised(2,:));
% title('Generated Random Pose Graph with Noise');
% xlabel('x'); ylabel('y');





