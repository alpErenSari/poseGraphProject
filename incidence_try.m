N = 10;

A = zeros(N,N);

% for i = 1:N
%     A(i, i+1) = 1;
%     A(i+1, i) = -1;
% end

for i = 1:N
    for j = 1:N
        prob = rand(1);
        if(prob < 0.1)
            A(i,j) = 1;
            A(j,i) = -1;
        end
    end
    if(i < N)
        A(i, i+1) = 1;
        A(i+1, i) = -1;
    end
end

pmax = 10;
p = pmax*rand(2,N);
theta = 2*pi*randn(1,N) - pi;
R_cells = euler_to_rot_mat(theta);
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

D_cell = cell(

    