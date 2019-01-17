function delta_ij = delta_p_cell_to_p_mat(p);
    delta_ij = cell(1,length(p));
    for i = 1:length(p)
        delta_x = p{i}(1);
        delta_y = p{i}(2);
        delta_ij{i} = [delta_x -delta_y; delta_y delta_x];
    end
end