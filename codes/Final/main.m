% 主程序
% 提取节点和单元信息
xn = msh.POS(:, 1); % 提取所有节点的 x 坐标
yn = msh.POS(:, 2); % 提取所有节点的 y 坐标
eln = msh.TRIANGLES(:, 1:3); % 提取所有单元的节点索引

% 材料属性
E = 1e10;  % 弹性模量 (Pa)
nu = 0.3;  % 泊松比
kappa = 1.0; % conductivity
% 平面应力本构矩阵
cm = (E / (1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2]; % 计算平面应力条件下的本构矩阵

% 自由度编号
nn = length(xn); % 总节点数
dof = reshape(1:(2 * nn), nn, 2); % 每个节点有两个自由度 (u, v)，重新排列为 [1, 2; 3, 4; ...]
ndof = 2 * nn; % 总方程数（总自由度数）

% LM 数组 (局部到全局自由度映射)
nelem = size(eln, 1); % 总单元数
npe = size(eln, 2); % 每个单元的节点数（3 个节点）
lm = zeros(nelem, 2 * npe); % 初始化 LM 数组，每个单元有 6 个自由度（3 节点 * 2 DOFs/节点）

% 向量化操作
for elem = 1:nelem % 遍历所有单元
    for node = 1:npe % 遍历每个单元的节点
        lm(elem, (node-1)*2 + 1:node*2) = dof(eln(elem, node), :); % 将局部自由度映射到全局自由度
    end
end

% 高斯积分点
gp = 2; % 每个方向的积分点数（2x2 高斯积分点）
[xi, eta, w] = Gauss2D(gp, gp); % 生成高斯积分点和权重

% 全局刚度矩阵和力向量
K = sparse(ndof, ndof); % 初始化全局刚度矩阵为稀疏矩阵
F = zeros(ndof, 1); % 初始化力向量

% 单元刚度矩阵组装
for elem = 1:nelem % 遍历所有单元
    xcoords = xn(eln(elem, :)); % 当前单元的 x 坐标
    ycoords = yn(eln(elem, :)); % 当前单元的 y 坐标
    k_elem = zeros(6, 6); % 初始化单元刚度矩阵
    for quad = 1:gp^2 % 遍历所有高斯积分点
        % 形函数及其导数
        N = zeros(1, 3); % 形函数矩阵
        Nx = zeros(1, 3); % 形函数对 xi 的导数
        Ny = zeros(1, 3); % 形函数对 eta 的导数
        for i = 1:3 % 遍历每个节点
            N(i) = Quad_tri(i, xi(quad), eta(quad)); % 计算形函数
            [Nx(i), Ny(i)] = Quad_grad_tri(i, xi(quad), eta(quad)); % 计算形函数的导数
        end
        % 映射
        J = [Nx * xcoords, Nx * ycoords; Ny * xcoords, Ny * ycoords]; % 计算雅可比矩阵
        [detJ, invJ] = detinv(J); % 计算雅可比矩阵的行列式和逆矩阵
        dNx = invJ(1, 1) * Nx + invJ(1, 2) * Ny; % 映射后的形函数对 x 的导数
        dNy = invJ(2, 1) * Nx + invJ(2, 2) * Ny; % 映射后的形函数对 y 的导数
        B = zeros(3, 6); % 初始化应变-位移矩阵
        for i = 1:3 % 遍历每个节点
            B(1, (i-1)*2 + 1) = dNx(i); % 应变-位移矩阵的第 1 行
            B(2, (i-1)*2 + 2) = dNy(i); % 应变-位移矩阵的第 2 行
            B(3, (i-1)*2 + 1:((i-1)*2 + 2)) = [dNy(i), dNx(i)]; % 应变-位移矩阵的第 3 行
        end
        % 计算单元刚度矩阵
        k_elem = k_elem + B' * cm * B * detJ * w(quad); % 积分求解单元刚度矩阵
    end
    % 组装到全局刚度矩阵
    for i = 1:6
        for j = 1:6
            if lm(elem, i) > 0 && lm(elem, j) > 0
                K(lm(elem, i), lm(elem, j)) = K(lm(elem, i), lm(elem, j)) + k_elem(i, j); % 将单元刚度矩阵组装到全局刚度矩阵
            end
        end
    end
end

% 边界条件
% 左边界 u = 0
left_nodes = find(abs(xn + 1) < 1e-5); % 找到左边界节点
for i = 1:length(left_nodes)
    K(dof(left_nodes(i), 1), :) = 0; % 在左边界节点的 u 方向设置刚度矩阵为 0
    K(dof(left_nodes(i), 1), dof(left_nodes(i), 1)) = 1; % 对角线设为 1
    F(dof(left_nodes(i), 1)) = 0; % 在左边界节点的 u 方向设置力向量为 0
end

% 下边界 v = 0
bottom_nodes = find(abs(yn + 1) < 1e-5); % 找到下边界节点
for i = 1:length(bottom_nodes)
    K(dof(bottom_nodes(i), 2), :) = 0; % 在下边界节点的 v 方向设置刚度矩阵为 0
    K(dof(bottom_nodes(i), 2), dof(bottom_nodes(i), 2)) = 1; % 对角线设为 1
    F(dof(bottom_nodes(i), 2)) = 0; % 在下边界节点的 v 方向设置力向量为 0
end

% 右边界施加拉伸力
right_nodes = find(abs(xn - 1) < 1e-5); % 找到右边界节点
force_per_node = 10e3 * 0.1; % 每个节点的拉伸力，10 KPa
for i = 1:length(right_nodes)
    F(dof(right_nodes(i), 1)) = F(dof(right_nodes(i), 1)) + force_per_node; % 在右边界节点的 u 方向施加拉伸力
end

% 求解位移
u = K \ F; % 使用稀疏矩阵求解线性方程组 K * u = F

% 位移后处理
u_dis = u(1:2:end); % 提取 u 方向的位移
v_dis = u(2:2:end); % 提取 v 方向的位移

% 计算单元应力和应变
stress = zeros(nelem, 3); % 初始化应力矩阵
strain = zeros(nelem, 3); % 初始化应变矩阵

for elem = 1:nelem % 遍历所有单元
    nodes = eln(elem, :); % 当前单元的节点索引
    xcoords = xn(nodes); % 当前单元的 x 坐标
    ycoords = yn(nodes); % 当前单元的 y 坐标
    d_elem = [u_dis(nodes); v_dis(nodes)]; % 当前单元的位移向量，维度 [6x1]
    elem_stress = zeros(3, 1); % 初始化当前单元的应力向量
    elem_strain = zeros(3, 1); % 初始化当前单元的应变向量
    
    for quad = 1:gp^2 % 遍历所有高斯积分点
        % 形函数及其导数
        N = zeros(1, 3); % 形函数矩阵
        Nx = zeros(1, 3); % 形函数对 xi 的导数
        Ny = zeros(1, 3); % 形函数对 eta 的导数
        for i = 1:3 % 遍历每个节点
            N(i) = Quad_tri(i, xi(quad), eta(quad)); % 计算形函数
            [Nx(i), Ny(i)] = Quad_grad_tri(i, xi(quad), eta(quad)); % 计算形函数的导数
        end
        % 映射
        J = [Nx * xcoords, Nx * ycoords; Ny * xcoords, Ny * ycoords]; % 计算雅可比矩阵
        [detJ, invJ] = detinv(J); % 计算雅可比矩阵的行列式和逆矩阵
        dNx = invJ(1, 1) * Nx + invJ(1, 2) * Ny; % 映射后的形函数对 x 的导数
        dNy = invJ(2, 1) * Nx + invJ(2, 2) * Ny; % 映射后的形函数对 y 的导数
        B = zeros(3, 6); % 初始化应变-位移矩阵
        for i = 1:3 % 遍历每个节点
            B(1, (i-1)*2 + 1) = dNx(i); % 应变-位移矩阵的第 1 行
            B(2, (i-1)*2 + 2) = dNy(i); % 应变-位移矩阵的第 2 行
            B(3, (i-1)*2 + 1:((i-1)*2 + 2)) = [dNy(i), dNx(i)]; % 应变-位移矩阵的第 3 行
        end
        % 计算高斯点应变
        strain_gauss = B * d_elem; % 计算高斯点的应变
        % 计算高斯点应力
        stress_gauss = cm * strain_gauss; % 计算高斯点的应力
        elem_stress = elem_stress + stress_gauss * w(quad); % 积分求解单元应力
        elem_strain = elem_strain + strain_gauss * w(quad); % 积分求解单元应变
    end
    % 平均应力和应变
    stress(elem, :) = elem_stress / sum(w); % 计算单元的平均应力
    strain(elem, :) = elem_strain / sum(w); % 计算单元的平均应变
end

% 计算等效应力 (Von Mises 应力)
von_mises = sqrt(stress(:, 1).^2 + stress(:, 2).^2 - stress(:, 1) .* stress(:, 2) + 3 * stress(:, 3).^2); % 计算 Von Mises 应力

% 可视化
% 应力
figure;
patch('Faces', eln, 'Vertices', [xn, yn], 'FaceVertexCData', von_mises, ... % 绘制单元应力分布
    'FaceColor', 'flat', 'EdgeColor', 'k');
colorbar;
title('四分之一带孔平板的应力分布 (Von Mises 应力)');
xlabel('X 坐标');
ylabel('Y 坐标');
axis equal;

% 应变
ex = strain(:, 1); % 提取 εx 应变
figure;
patch('Faces', eln, 'Vertices', [xn, yn], 'FaceVertexCData', ex, ... % 绘制单元 εx 应变分布
    'FaceColor', 'flat', 'EdgeColor', 'k');
colorbar;
title('四分之一带孔平板的应变分布 (εx)');
xlabel('X 坐标');
ylabel('Y 坐标');
axis equal;

% 辅助函数
function [detJ, invJ] = detinv(J)
    detJ = det(J);
    invJ = inv(J);
end

function val = Quad_tri(aa, xi, eta)
    if aa == 1
        val = 1 - xi - eta;
    elseif aa == 2
        val = xi;
    elseif aa == 3
        val = eta;
    else
        error('Error: value of a should be 1, 2, or 3.');
    end
end

function [Na_xi, Na_eta] = Quad_grad_tri(aa, xi, eta)
    if aa == 1
        Na_xi = -1;
        Na_eta = -1;
    elseif aa == 2
        Na_xi = 1;
        Na_eta = 0;
    elseif aa == 3
        Na_xi = 0;
        Na_eta = 1;
    else
        error('Error: value of a should be 1, 2, or 3.');
    end
end

function [xi, eta, w] = Gauss2D(N1, N2)
    % 预分配
    xi = zeros(N1 * N2, 1);
    eta = xi;
    w = xi;
    
    % 生成1D规则
    [x1, w1] = Gauss(N1, -1, 1);
    [x2, w2] = Gauss(N2, -1, 1);
    
    k = 1;
    for i = 1:N1
        for j = 1:N2
            xi(k) = x1(i);
            eta(k) = x2(j);
            w(k) = w1(i) * w2(j);
            k = k + 1;
        end
    end
end

% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
function [x,w] = Gauss(N,a,b)
N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;   
     
% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end
