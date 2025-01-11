function disp = solve_system(K, F, nodes, elements, boundary_edges)
    num_nodes = size(nodes, 1);
    num_dofs = 2 * num_nodes;

    % 处理位移边界条件
    fixed_dofs = [];
    for be = 1:size(boundary_edges, 1)
        n1 = boundary_edges(be, 1);
        n2 = boundary_edges(be, 2);
        x1 = nodes(n1, 1);
        y1 = nodes(n1, 2);
        x2 = nodes(n2, 1);
        y2 = nodes(n2, 2);
        if y1 == 0 && y2 == 0 % 底部边界（对称边界）
            fixed_dofs = [fixed_dofs; 2 * n1 - 1; 2 * n2 - 1];
        elseif x1 == 0 && x2 == 0 % 左边界（固定边界）
            fixed_dofs = [fixed_dofs; 2 * n1 - 1; 2 * n1; 2 * n2 - 1; 2 * n2];
        end
    end

    % 从刚度矩阵和载荷向量中删除固定自由度
    free_dofs = setdiff(1:num_dofs, fixed_dofs);
    K_free = K(free_dofs, free_dofs);
    F_free = F(free_dofs);

    % 求解
    disp_free = K_free \ F_free;
    disp = zeros(num_dofs, 1);
    disp(free_dofs) = disp_free;
end