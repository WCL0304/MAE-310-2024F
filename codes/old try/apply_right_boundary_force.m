% 右边界施加拉伸力
function apply_right_boundary_force(F, x_coor, ID)
    right_boundary_nodes = find(abs(x_coor - 1) < 1e-5);
    force_per_node = 10e3 * 0.1; % 10 KPa
    for i = 1:length(right_boundary_nodes)
        dof = ID(right_boundary_nodes(i), 1);
        F(dof) = F(dof) + force_per_node;
    end
end