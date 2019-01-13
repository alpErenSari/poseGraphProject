function cost = cost_calculator(p, R, p_delta, R_deltas)
    r_vec = rot_mat_to_vec(R);
    delta_p_mat = delta_p_to_p_mat(p_delta);
    cost = 0;
    for i = 1:length(p)-1
        cost = cost + norm(p{i+1}-p{i} - delta_p_mat{i}*r_vec{i}).^2 ...
            + norm(r_vec{i+1} - R_deltas{i}*r_vec{i}).^2;
    end
end