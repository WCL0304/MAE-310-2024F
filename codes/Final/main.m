% 主程序
% 提取节点和单元信息
x = msh.POS(:, 1); % 节点的 x 坐标
y = msh.POS(:, 2); % 节点的 y 坐标
elem_nodes = msh.TRIANGLES(:, 1:3); % 单元的节点索引
% 材料属性
E = 1e10;  % 弹性模量 (Pa)
nu = 0.3;  % 泊松比
% 平面应力本构矩阵
D = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
% 自由度编号
num_nodes = length(x); % 总节点数
dof = reshape(1:(2 * num_nodes), num_nodes, 2); % 每个节点两个自由度 (u, v)
total_dof = 2 * num_nodes; % 总自由度数
% LM 数组 (局部到全局自由度映射)
num_elem = size(elem_nodes, 1); % 总单元数
nodes_per_elem = 3; % 每个单元的节点数
lm = zeros(num_elem, 2 * nodes_per_elem); % 单元自由度映射数组
% 向量化操作
for node_idx = 1:nodes_per_elem
    lm(:, (node_idx-1)*2 + 1:node_idx*2) = dof(elem_nodes(:, node_idx), :);
end
% 高斯积分点
gp = 2; % 每个方向的积分点数
[xi, eta, w] = Gauss2D(gp, gp);
% 全局刚度矩阵和力向量
K = sparse(total_dof, total_dof); % 初始化全局刚度矩阵
F = zeros(total_dof, 1); % 初始化力向量