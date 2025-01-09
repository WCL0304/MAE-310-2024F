% 主程序
clear; clc;

% 读取网格数据
[nodes, elements, boundaryEdges] = read_gmsh_mesh('quarter-plate-with-hole-quad.msh');
disp(['节点数组大小:', num2str(size(nodes, 1)), 'x', num2str(size(nodes, 2))]); % 调试信息
disp(['边界元素数组大小:', num2str(size(boundaryEdges, 1)), 'x', num2str(size(boundaryEdges, 2))]); % 调试信息
disp('边界元素内容:');
disp(boundaryEdges); % 打印边界元素的具体内容

% 材料参数
E = 1e11; % 弹性模量
nu = 0.3; % 泊松比
L = 4;    % 假设L为4
R = 0.5;  % 孔的半径

% 模型选择
isPlaneStress = true; % 平面应力模型

% 处理边界条件
[dirichletBC, neumannBC] = apply_boundary_conditions(boundaryEdges, nodes, L);

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

% 应用Neumann边界条件
for i = 1:size(neumannBC, 1)
    node = neumannBC(i, 1);
    direction = neumannBC(i, 2);
    value = neumannBC(i, 3);
    F_global(2 * (node - 1) + direction) = value;
end

% 应用Dirichlet边界条件
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
exact_solution = @(x, y, type) exact_solution_helper(x, y, L, R, type);
[L2_error, H1_error] = calculate_errors(nodes, displacements, exact_solution);
disp(['L2范数误差:', num2str(L2_error)]);
disp(['H1范数误差:', num2str(H1_error)]);