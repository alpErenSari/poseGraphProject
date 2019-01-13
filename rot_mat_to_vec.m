function rot_vec = rot_mat_to_vec(R)
    rot_vec = cell(1,length(R));
    for i = 1:length(R)
        rot_vec{i} = R{i}(:,1);
    end
end