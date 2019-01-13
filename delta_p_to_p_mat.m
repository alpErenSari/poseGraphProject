function delta_ij = delta_p_to_p_mat(p);
    delta_ij = cell(1,length(p));
    for i = 1:length(p)
        delta_x = p(1,i);
        delta_y = p(2,i);
        delta_ij{i} = [delta_x -delta_y; delta_y delta_x];
    end
end