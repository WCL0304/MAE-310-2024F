function applyBoundaryConditions(x, y, dof, K, F)
    % 施加边界条件

    % 左边界 u = 0
    left_boundary_nodes = find(abs(x + 1) < 1e-5);
    num_left_nodes = numel(left_boundary_nodes);
    for i = 1:num_left_nodes
        u_dof = dof(left_boundary_nodes(i), 1);
        K(u_dof, :) = 0;
        K(u_dof, u_dof) = 1;
        F(u_dof) = 0;
    end

    % 下边界 v = 0
    bottom_boundary_nodes = find(abs(y + 1) < 1e-5);
    num_bottom_nodes = numel(bottom_boundary_nodes);
    for i = 1:num_bottom_nodes
        v_dof = dof(bottom_boundary_nodes(i), 2);
        K(v_dof, :) = 0;
        K(v_dof, v_dof) = 1;
        F(v_dof) = 0;
    end

    % 右边界施加拉伸力
    right_boundary_nodes = find(abs(x - 1) < 1e-5);
    num_right_nodes = numel(right_boundary_nodes);
    force_per_node = 10e3 * 0.1; % 10 KPa
    for i = 1:num_right_nodes
        u_dof = dof(right_boundary_nodes(i), 1);
        F(u_dof) = F(u_dof) + force_per_node;
    end
end