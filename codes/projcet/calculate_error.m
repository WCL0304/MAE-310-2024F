function [error_L2, error_H1] = calculate_error(disp, nodes, elements, exact_disp, exact_strain, manufact_solution, manufact_solution_grad_x, manufact_solution_grad_y)
    num_nodes = size(nodes, 1);
    num_elements = size(elements, 1);
    disp_x = disp(1:2:end);
    disp_y = disp(2:2:end);

    u_num = [disp_x, disp_y];
    u_exact = exact_disp(nodes(:, 1), nodes(:, 2));

    % 位移误差
    error_L2 = sqrt(sum((u_num - u_exact).^2) / sum(u_exact.^2));

    % 应变误差
    strain_num = zeros(num_nodes, 3);
    for ee = 1:num_elements
        x_ele = nodes(elements(ee, :), 1);
        y_ele = nodes(elements(ee, :), 2);
        disp_e = [disp_x(elements(ee, :)); disp_y(elements(ee, :))];
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

            strain = B * disp_e;
            strain_num(elements(ee, :), :) = strain_num(elements(ee, :), :) + strain * weight(ll) * detJ;
        end
    end

    strain_exact = exact_strain(nodes(:, 1), nodes(:, 2));
    error_H1 = sqrt(sum((strain_num - strain_exact).^2) / sum(strain_exact.^2));
end