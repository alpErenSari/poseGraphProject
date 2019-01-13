function rot_mat = euler_to_rot_mat(theta)
%     rot_mat = zeros(2,2*length(theta));
    rot_mat = cell(1,length(theta));
    for i = 1:length(theta)
        rot_mat{i} = [cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))];
    end
end