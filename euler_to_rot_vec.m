function rot_vec = euler_to_rot_vec(theta)
rot_vec = cell(1,length(theta));
    for i = 1:length(theta)
        rot_vec{i} = [cos(theta(i)) sin(theta(i))]';
    end
end
    