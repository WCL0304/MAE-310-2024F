% 单元刚度矩阵组装
function K = assemble_global_stiffness_matrix(K, IEN, x_coor, y_coor, C, xi, eta, weight)
    n_elements = size(IEN, 1);
    for ee = 1:n_elements
        x_nodes = x_coor(IEN(ee, :));
        y_nodes = y_coor(IEN(ee, :));
        k_elem = compute_element_stiffness_matrix(x_nodes, y_nodes, C, xi, eta, weight);
        K = add_element_stiffness_to_global(K, IEN(ee, :), k_elem);
    end
end
