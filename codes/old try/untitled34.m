% 主程序
% 提取节点和单元信息
x_coor = msh.POS(:, 1); % x 坐标
y_coor = msh.POS(:, 2); % y 坐标
IEN = msh.TRIANGLES(:, 1:3); % 三角形单元 (节点索引)

% 材料参数
E = 1e10; % 弹性模量（Pa）
nu = 0.3; % 泊松比
L = 4; % 假设 L 为 4

% 模型选择
isPlaneStress = true; % 平面应力模型

% 处理边界条件
dirichletBC = [];
neumannBC = [];
left_boundary_nodes = find(abs(nodes(:, 1) - (-1)) < 1e-5);
dirichletBC = [dirichletBC; left_boundary_nodes, 1, 0]; % 左边界 u = 0

bottom_boundary_nodes = find(abs(nodes(:, 2) - (-1)) < 1e-5);
dirichletBC = [dirichletBC; bottom_boundary_nodes, 2, 0]; % 下边界 v = 0

right_boundary_nodes = find(abs(nodes(:, 1) - 1) < 1e-5);
T = 10e3 * 0.1; % 10 KPa
neumannBC = [neumannBC; right_boundary_nodes, 1, T]; % 右边界施加拉伸力

% 初始化全局刚度矩阵和载荷向量
K_global = sparse(2 * size(nodes, 1), 2 * size(nodes, 1));
F_global = zeros(2 * size(nodes, 1), 1);

% 集成刚度矩阵和载荷向量
for i = 1:size(elements, 1)
    element = elements(i, :);
    K_elem = calculate_element_stiffness_matrix(element, nodes, E, nu, isPlaneStress);
    indices = 2 * (element - 1) + [1, 2];
    K_global(indices, indices) = K_global(indices, indices) + K_elem;
end

% 应用 Neumann 边界条件
for i = 1:size(neumannBC, 1)
    node = neumannBC(i, 1);
    direction = neumannBC(i, 2);
    value = neumannBC(i, 3);
    F_global(2 * (node - 1) + direction) = value;
end

% 应用 Dirichlet 边界条件
K_dirichlet = K_global;
F_dirichlet = F_global;
for i = 1:size(dirichletBC, 1)
    node = dirichletBC(i, 1);
    direction = dirichletBC(i, 2);
    value = dirichletBC(i, 3);
    K_dirichlet(2 * (node - 1) + direction, :) = 0;
    K_dirichlet(:, 2 * (node - 1) + direction) = 0;
    K_dirichlet(2 * (node - 1) + direction, 2 * (node - 1) + direction) = 1;
    F_dirichlet(2 * (node - 1) + direction) = value;
end

% 求解位移
displacements = K_dirichlet \ F_dirichlet;

% 计算应变和应力
strains = calculate_strains(displacements, nodes, elements);
stresses = calculate_stresses(strains, E, nu, isPlaneStress);

% 可视化结果
visualize_results(nodes, displacements, strains, stresses);

% 计算误差
exact_solution = @(x, y, type) exact_solution_helper(x, y, L, 1, type);
[L2_error, H1_error] = calculate_errors(nodes, displacements, exact_solution);
disp(['L2范数误差:', num2str(L2_error)]);
disp(['H1范数误差:', num2str(H1_error)]);

% 计算单元刚度矩阵
function K_elem = calculate_element_stiffness_matrix(element, nodes, E, nu, isPlaneStress)
    % 计算平面应力/应变矩阵
    D = plane_stress_strain_matrix(E, nu, isPlaneStress);
    
    % 高斯积分点
    n_int = 2; % 每个方向的积分点数
    [xi, eta, weight] = Gauss2D(n_int, n_int);
    
    % 获取单元节点坐标
    x_nodes = nodes(element, 1);
    y_nodes = nodes(element, 2);
    
    k_elem = zeros(6, 6); % 单元刚度矩阵
    for ll = 1:n_int^2
        % 形函数及其导数
        N = zeros(1, 3);
        dN_dxi = zeros(1, 3);
        dN_deta = zeros(1, 3);
        for aa = 1:3
            N(aa) = Quad_tri(aa, xi(ll), eta(ll));
            [dN_dxi(aa), dN_deta(aa)] = Quad_grad_tri(aa, xi(ll), eta(ll));
        end
        
        % 映射到物理域
        J = [dN_dxi * x_nodes, dN_dxi * y_nodes; dN_deta * x_nodes, dN_deta * y_nodes];
        detJ = det(J);
        invJ = inv(J);
        dN_dx = invJ(1, 1) * dN_dxi + invJ(1, 2) * dN_deta;
        dN_dy = invJ(2, 1) * dN_dxi + invJ(2, 2) * dN_deta;
        
        B = zeros(3, 6); % 应变-位移矩阵
        for i = 1:3
            B(:, (i-1)*2 + 1) = [dN_dx(i); 0; dN_dy(i)];
            B(:, (i-1)*2 + 2) = [0; dN_dy(i); dN_dx(i)];
        end
        
        % 计算单元刚度矩阵
        k_elem = k_elem + B' * D * B * detJ * weight(ll);
    end
    K_elem = k_elem;
end

% 计算应变
function strains = calculate_strains(displacements, nodes, elements)
    nElements = size(elements, 1);
    strains = zeros(nElements, 3);
    for i = 1:nElements
        element = elements(i, :);
        disp_elem = displacements(2 * (element - 1) + [1, 2], :);
        B = strain_displacement_matrix(element, nodes, [0, 0]);
        strains(i, :) = B * disp_elem / det(B);
    end
end

% 计算应力
function stresses = calculate_stresses(strains, E, nu, isPlaneStress)
    D = plane_stress_strain_matrix(E, nu, isPlaneStress);
    stresses = D * strains';
end

% 计算平面应力/应变矩阵
function D = plane_stress_strain_matrix(E, nu, isPlaneStress)
    if isPlaneStress
        D = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    else
        D = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
    end
end

% 计算应变-位移矩阵和雅各比矩阵的行列式
function B = strain_displacement_matrix(element, nodes, gauss_point)
    % 假设使用线性三角形单元
    p1 = nodes(element(1), :);
    p2 = nodes(element(2), :);
    p3 = nodes(element(3), :);
    
    % 形函数及其导数
    N = zeros(1, 3);
    dN_dxi = zeros(1, 3);
    dN_deta = zeros(1, 3);
    for aa = 1:3
        N(aa) = Quad_tri(aa, gauss_point(1), gauss_point(2));
        [dN_dxi(aa), dN_deta(aa)] = Quad_grad_tri(aa, gauss_point(1), gauss_point(2));
    end
    
    % 雅各比矩阵
    J = [dN_dxi * [p1(1); p2(1); p3(1)], dN_dxi * [p1(2); p2(2); p3(2)]; ...
        dN_deta * [p1(1); p2(1); p3(1)], dN_deta * [p1(2); p2(2); p3(2)]];
    detJ = det(J);
    invJ = inv(J);
    dN_dx = invJ(1, 1) * dN_dxi + invJ(1, 2) * dN_deta;
    dN_dy = invJ(2, 1) * dN_dxi + invJ(2, 2) * dN_deta;
    
    % 应变-位移矩阵
    B = zeros(3, 6);
    for i = 1:3
        B(:, (i-1)*2 + 1) = [dN_dx(i); 0; dN_dy(i)];
        B(:, (i-1)*2 + 2) = [0; dN_dy(i); dN_dx(i)];
    end
end

% 可视化结果
function visualize_results(nodes, displacements, strains, stresses)
    % 绘制位移
    figure;
    quiver(nodes(:, 1), nodes(:, 2), displacements(1:2:end), displacements(2:2:end));
    title('位移分布图');
    xlabel('x');
    ylabel('y');
    axis equal;
    
    % 绘制应变
    figure;
    contourf(nodes(:, 1), nodes(:, 2), strains);
    colorbar;
    title('应变分布图');
    xlabel('x');
    ylabel('y');
    axis equal;
    
    % 绘制应力
    figure;
    contourf(nodes(:, 1), nodes(:, 2), stresses);
    colorbar;
    title('应力分布图');
    xlabel('x');
    ylabel('y');
    axis equal;
end

% 计算解析解
function [u, v] = exact_solution_helper(x, y, L, R, type)
    T = 10e3; % 牵引力
    r = sqrt(x.^2 + y.^2);
    theta = atan2(y, x);
    if strcmp(type, 'u')
        u = T * (1 - R^2 ./ r^2) * cos(theta);
        v = 0; % 只计算 u 方向的位移
    elseif strcmp(type, 'v')
        u = 0; % 只计算 v 方向的位移
        v = T * (1 - R^2 ./ r^2) * sin(theta);
    end
end

% 计算误差
function [L2_error, H1_error] = calculate_errors(nodes, displacements, exact_solution)
    nNodes = size(nodes, 1);
    x_coor = nodes(:, 1);
    y_coor = nodes(:, 2);
    u = displacements(1:2:end);
    v = displacements(2:2:end);
    
    % 计算解析解
    u_exact = zeros(nNodes, 1);
    v_exact = zeros(nNodes, 1);
    for i = 1:nNodes
        [u_exact(i), ~] = exact_solution(x_coor(i), y_coor(i), 'u');
        [~, v_exact(i)] = exact_solution(x_coor(i), y_coor(i), 'v');
    end
    
    % L2 范数误差
    L2_error_u = norm(u - u_exact) / norm(u_exact);
    L2_error_v = norm(v - v_exact) / norm(v_exact);
    L2_error = mean([L2_error_u, L2_error_v]);
    
    % 计算梯度
    [dudx, dudy] = gradient(u, x_coor, y_coor);
    [dvdx, dvdy] = gradient(v, x_coor, y_coor);
    [dudx_exact, dudy_exact] = gradient(u_exact, x_coor, y_coor);
    [dvdx_exact, dvdy_exact] = gradient(v_exact, x_coor, y_coor);
    
    % H1 范数误差
    H1_error_u = norm([dudx - dudx_exact; dudy - dudy_exact]) / norm([dudx_exact; dudy_exact]);
    H1_error_v = norm([dvdx - dvdx_exact; dvdy - dvdy_exact]) / norm([dvdx_exact; dvdy_exact]);
    H1_error = mean([H1_error_u, H1_error_v]);
end