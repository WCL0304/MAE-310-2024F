function [K, F] = assemble_matrix_vector(nodes, elements, traction, stress_type)
    num_nodes = size(nodes, 1);
    num_elements = size(elements, 1);
    K = sparse(2 * num_nodes, 2 * num_nodes);
    F = zeros(2 * num_nodes, 1);

    n_int_xi = 2;
    n_int_eta = 2;
    [xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

    for ee = 1:num_elements
        x_ele = nodes(elements(ee, :), 1);
        y_ele = nodes(elements(ee, :), 2);

        k_ele = zeros(2 * size(elements, 2), 2 * size(elements, 2));
        f_ele = zeros(2 * size(elements, 2), 1);

        for ll = 1:length(weight)
            if size(elements, 2) == 4
                [dN_dxi, dN_deta] = Quad_grad(1:4, xi(ll), eta(ll));
                N = Quad(1:4, xi(ll), eta(ll));
            elseif size(elements, 2) == 3
                [dN_dxi, dN_deta] = Tri_grad(1:3, xi(ll), eta(ll));
                N = Tri(1:3, xi(ll), eta(ll));
            end

            x = x_ele' * N;
            y = y_ele' * N;

            dx_dxi = dN_dxi * x_ele;
            dx_deta = dN_deta * x_ele;
            dy_dxi = dN_dxi * y_ele;
            dy_deta = dN_deta * y_ele;

            detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
            B = [dN_dxi(1) 0 dN_dxi(2) 0 dN_dxi(3) 0 dN_dxi(4) 0;
                  0 dN_deta(1) 0 dN_deta(2) 0 dN_deta(3) 0 dN_deta(4);
                  dN_deta(1) dN_dxi(1) dN_deta(2) dN_dxi(2) dN_deta(3) dN_dxi(3) dN_deta(4) dN_dxi(4)];

            if strcmp(stress_type, 'plane_stress')
                D = kappa * [1 nu 0; nu 1 0; 0 0 0];
            elseif strcmp(stress_type, 'plane_strain')
                D = 2 * mu * [1 -nu 0; -nu 1 0; 0 0 1 -nu] + lambda * [1 1 0; 1 1 0; 0 0 0];
            else
                error('未知的应力类型：%s', stress_type);
            end

            k_ele = k_ele + weight(ll) * detJ * (B' * D * B);
            f_ele = f_ele + weight(ll) * detJ * traction * [N'; N'];
        end

        % 将元素刚度矩阵和载荷向量组装到全局刚度矩阵和载荷向量
        for aa = 1:size(elements, 2)
            for bb = 1:size(elements, 2)
                global_dof_aa = 2 * elements(ee, aa) - 1:2 * elements(ee, aa);
                global_dof_bb = 2 * elements(ee, bb) - 1:2 * elements(ee, bb);
                K(global_dof_aa, global_dof_bb) = K(global_dof_aa, global_dof_bb) + k_ele(2 * (aa - 1) + 1:2 * aa, 2 * (bb - 1) + 1:2 * bb);
            end
        end

        global_dof = 2 * elements(ee, :) - 1:2 * elements(ee, :);
        F(global_dof) = F(global_dof) + f_ele;
    end

    % 处理边界条件
    for be = 1:size(boundary_edges, 1)
        n1 = boundary_edges(be, 1);
        n2 = boundary_edges(be, 2);
        x1 = nodes(n1, 1);
        y1 = nodes(n1, 2);
        x2 = nodes(n2, 1);
        y2 = nodes(n2, 2);
        length = norm([x2 - x1; y2 - y1]);

        if y1 == 0 && y2 == 0 % 底部边界（对称边界）
            F(2 * n1 - 1) = 0; % 对称边界，x 方向位移为 0
            F(2 * n2 - 1) = 0;
        elseif x1 == 4 && x2 == 4 % 右边界（Neumann 边界）
            F(2 * n1) = F(2 * n1) + traction * length * 0.5;
            F(2 * n2) = F(2 * n2) + traction * length * 0.5;
        end
    end

    % 处理固定边界条件（对称边界）
    fixed_dofs = [];
    for be = 1:size(boundary_edges, 1)
        n1 = boundary_edges(be, 1);
        n2 = boundary_edges(be, 2);
        x1 = nodes(n1, 1);
        y1 = nodes(n1, 2);
        x2 = nodes(n2, 1);
        y2 = nodes(n2, 2);
        if y1 == 0 && y2 == 0 % 底部边界
            fixed_dofs = [fixed_dofs; 2 * n1 - 1; 2 * n2 - 1];
        elseif x1 == 0 && x2 == 0 % 左边界
            fixed_dofs = [fixed_dofs; 2 * n1 - 1; 2 * n1; 2 * n2 - 1; 2 * n2];
        end
    end

    % 从刚度矩阵和载荷向量中删除固定自由度
    free_dofs = setdiff(1:2 * num_nodes, fixed_dofs);
    K_free = K(free_dofs, free_dofs);
    F_free = F(free_dofs);

    % 求解
    disp_free = K_free \ F_free;
    disp = zeros(2 * num_nodes, 1);
    disp(free_dofs) = disp_free;
end