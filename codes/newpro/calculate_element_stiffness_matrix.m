function K = calculate_element_stiffness_matrix(element, nodes, E, nu, isPlaneStress)
    % 计算单元刚度矩阵
    nNodes = size(element, 2);
    D = plane_stress_strain_matrix(E, nu, isPlaneStress);
    B = strain_displacement_matrix(element, nodes);

    if nNodes == 3 % 三角形单元
        area = polyarea(nodes(element, 1), nodes(element, 2));
        K = area * B' * D * B;
    elseif nNodes == 4 % 四边形单元
        gaussPoints = [-0.57735, -0.57735; 0.57735, -0.57735; 0.57735, 0.57735; -0.57735, 0.57735];
        K = zeros(2*nNodes, 2*nNodes);
        for i = 1:size(gaussPoints, 1)
            [J, detJ] = jacobian(element, nodes, gaussPoints(i, :));
            B = strain_displacement_matrix(element, nodes, gaussPoints(i, :));
            K = K + detJ * B' * D * B;
        end
    end
end

function D = plane_stress_strain_matrix(E, nu, isPlaneStress)
    % 计算弹性矩阵 D
    if isPlaneStress
        D = (E / (1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    else
        D = (E / ((1 + nu) * (1 - 2 * nu))) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
    end
end

function [B, detJ] = strain_displacement_matrix(element, nodes, xi)
    % 计算应变-位移矩阵 B 和雅克比行列式 detJ
    nNodes = size(element, 2);
    N = shape_functions(xi, nNodes);
    dNdx = shape_function_gradients(xi, nNodes, nodes(element, :));
    B = zeros(3, 2*nNodes);
    for i = 1:nNodes
        B(1, 2*i-1) = dNdx(i, 1);
        B(2, 2*i) = dNdx(i, 2);
        B(3, 2*i-1) = dNdx(i, 2);
        B(3, 2*i) = dNdx(i, 1);
    end

    [J, detJ] = jacobian(element, nodes, xi);
    B = B / detJ;
end

function [J, detJ] = jacobian(element, nodes, xi)
    % 计算雅克比矩阵 J 和其行列式 detJ
    nNodes = size(element, 2);
    N = shape_functions(xi, nNodes);
    dNdr = shape_function_gradients(xi, nNodes, nodes(element, :));
    J = [nodes(element, 1) * dNdr(:, 1), nodes(element, 2) * dNdr(:, 1); ...
         nodes(element, 1) * dNdr(:, 2), nodes(element, 2) * dNdr(:, 2)];
    detJ = det(J);
end

function N = shape_functions(xi, nNodes)
    % 计算形状函数 N
    if nNodes == 3 % 三角形单元
        N = [1 - xi(1) - xi(2); xi(1); xi(2)];
    elseif nNodes == 4 % 四边形单元
        N = [1/4 * (1 - xi(1)) * (1 - xi(2)); 1/4 * (1 + xi(1)) * (1 - xi(2)); ...
             1/4 * (1 + xi(1)) * (1 + xi(2)); 1/4 * (1 - xi(1)) * (1 + xi(2))];
    end
end

function dNdr = shape_function_gradients(xi, nNodes, nodes)
    % 计算形状函数梯度 dNdr
    if nNodes == 3 % 三角形单元
        dNdr = [-1; 1; 0; -1; 0; 1];
    elseif nNodes == 4 % 四边形单元
        dNdr = [1/4 * (1 - xi(2)), 1/4 * (1 - xi(1)); 1/4 * (1 + xi(2)), -1/4 * (1 + xi(1)); ...
                1/4 * (1 + xi(2)), 1/4 * (1 + xi(1)); -1/4 * (1 - xi(2)), 1/4 * (1 - xi(1))];
    end
end