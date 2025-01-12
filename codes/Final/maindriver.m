% 主程序
% 提取节点和单元信息
xn = msh.POS(:, 1); % 提取所有节点的 x 坐标
yn = msh.POS(:, 2); % 提取所有节点的 y 坐标
eln = msh.TRIANGLES(:, 1:3); % 提取所有单元的节点索引

% 材料属性
E = 1e9;  % 弹性模量 (Pa)
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
        N = arrayfun(@(i) Quad_tri(i, xi(quad), eta(quad)), 1:3); % 计算形函数
        [Nx, Ny] = deal(zeros(1, 3));
        for i = 1:3
            [Nx(i), Ny(i)] = Quad_grad_tri(i, xi(quad), eta(quad)); % 计算形函数的导数
        end
        
        % 映射
        J = [Nx * xcoords, Nx * ycoords; Ny * xcoords, Ny * ycoords]; % 计算雅可比矩阵
        [detJ, invJ] = detinv(J); % 计算雅可比矩阵的行列式和逆矩阵
        dNx = invJ(1, 1) * Nx + invJ(1, 2) * Ny; % 映射后的形函数对 x 的导数
        dNy = invJ(2, 1) * Nx + invJ(2, 2) * Ny; % 映射后的形函数对 y 的导数
        
        % 应变-位移矩阵
        B = zeros(3, 6);
        for i = 1:3
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
    idx = dof(left_nodes(i), 1);
    K(idx, :) = 0; % 在左边界节点的 u 方向设置刚度矩阵为 0
    K(idx, idx) = 1; % 对角线设为 1
    F(idx) = 0; % 在左边界节点的 u 方向设置力向量为 0
end

% 下边界 v = 0
bottom_nodes = find(abs(yn + 1) < 1e-5); % 找到下边界节点
for i = 1:length(bottom_nodes)
    idx = dof(bottom_nodes(i), 2);
    K(idx, :) = 0; % 在下边界节点的 v 方向设置刚度矩阵为 0
    K(idx, idx) = 1; % 对角线设为 1
    F(idx) = 0; % 在下边界节点的 v 方向设置力向量为 0
end

% 右边界施加拉伸力
right_nodes = find(abs(xn - 1) < 1e-5); % 找到右边界节点
force_per_node = 10e3 * 0.1; % 每个节点的拉伸力，10 KPa
for i = 1:length(right_nodes)
    idx = dof(right_nodes(i), 1);
    F(idx) = F(idx) + force_per_node; % 在右边界节点的 u 方向施加拉伸力
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
        N = arrayfun(@(i) Quad_tri(i, xi(quad), eta(quad)), 1:3); % 计算形函数
        [Nx, Ny] = deal(zeros(1, 3));
        for i = 1:3
            [Nx(i), Ny(i)] = Quad_grad_tri(i, xi(quad), eta(quad)); % 计算形函数的导数
        end
        
        % 映射
        J = [Nx * xcoords, Nx * ycoords; Ny * xcoords, Ny * ycoords]; % 计算雅可比矩阵
        [detJ, invJ] = detinv(J); % 计算雅可比矩阵的行列式和逆矩阵
        dNx = invJ(1, 1) * Nx + invJ(1, 2) * Ny; % 映射后的形函数对 x 的导数
        dNy = invJ(2, 1) * Nx + invJ(2, 2) * Ny; % 映射后的形函数对 y 的导数
        
        % 应变-位移矩阵
        B = zeros(3, 6);
        for i = 1:3
            B(1, (i-1)*2 + 1) = dNx(i); % 应变-位移矩阵的第 1 行
            B(2, (i-1)*2 + 2) = dNy(i); % 应变-位移矩阵的第 2 行
            B(3, (i-1)*2 + 1:((i-1)*2 + 2)) = [dNy(i), dNx(i)]; % 应变-位移矩阵的第 3 行
        end
        
        % 计算高斯点应变和应力
        strain_gauss = B * d_elem; % 计算高斯点的应变
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
patch('Faces', eln, 'Vertices', [xn, yn], 'FaceVertexCData', von_mises, ...
    'FaceColor', 'flat', 'EdgeColor', 'k');
colorbar;
title('四分之一带孔平板的应力分布 (Von Mises 应力)');
xlabel('X 坐标');
ylabel('Y 坐标');
axis equal;

% 应变
ex = strain(:, 1); % 提取 εx 应变
figure;
patch('Faces', eln, 'Vertices', [xn, yn], 'FaceVertexCData', ex, ...
    'FaceColor', 'flat', 'EdgeColor', 'k');
colorbar;
title('四分之一带孔平板的应变分布 (εx)');
xlabel('X 坐标');
ylabel('Y 坐标');
axis equal;

% % exact solution
exact = @(x,y) x.*(1-x).*y.*(1-y);
exact_x = @(x,y) (1-2*x).*y.*(1-y);
exact_y = @(x,y) x.*(1-x).*(1-2*y);

f = @(x,y) 2.0*kappa*x.*(1-x) + 2.0*kappa*y.*(1-y); % source term

% % quadrature rule
n_int_xi  = 2;
n_int_eta = 2;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% % mesh generation
n_el_x_values = msh.POS(:, 1);
n_el_y_values = msh.POS(:, 2);

% % 存储误差
errors_L2_quad = zeros(size(n_el_x_values));
errors_H1_quad = zeros(size(n_el_x_values));
errors_L2_tri = zeros(size(n_el_x_values));
errors_H1_tri = zeros(size(n_el_x_values));

% % 制造的解（不用注释）
manufactured_solution = @(x, y) x.^2 .* y.^2;

% % 制造的解的梯度（不用注释）
manufact_solution_grad_x = @(x, y) 2 * x .* y.^2;
manufact_solution_grad_y = @(x, y) 2 * y .* x.^2;

% % 计算不同网格尺寸下的误差（不用注释）
for i = 1:length(n_el_x_values)
    n_el_x = n_el_x_values(i);% 当前网格的 x 方向元素数
    n_el_y = n_el_y_values(i);% 当前网格的 y 方向元素数
    n_el = n_el_x * n_el_y;
    
    n_np_x = n_el_x + 1;% x 方向的节点数
    n_np_y = n_el_y + 1;% y 方向的节点数
    n_np = n_np_x * n_np_y;
    
    hx = 1.0 / n_el_x;
    hy = 1.0 / n_el_y;
    
    % % 生成节点坐标（不用注释）
    x_coor = zeros(n_np, 1);
    y_coor = x_coor;
    for ny = 1 : n_np_y
        for nx = 1 : n_np_x
            index = (ny-1)*n_np_x + nx;
            x_coor(index) = (nx-1) * hx;
            y_coor(index) = (ny-1) * hy;
        end
    end
    % 重新定义IEN数组和ID数组（三角形代码）
    IEN_tri = zeros(2*n_el, 3);
    for ex = 1 : n_el_x
        for ey = 1 : n_el_y
            ee = (ey-1) * n_el_x + ex;
            IEN_tri(2*(ee-1)+1, :) = [(ey-1) * n_np_x + ex, (ey-1) * n_np_x + ex + 1, ey * n_np_x + ex + 1]; % 第一个三角形元素
            IEN_tri(2*ee, :) = [(ey-1) * n_np_x + ex, ey * n_np_x + ex + 1, ey * n_np_x + ex]; % 第二个三角形元素
        end
    end

    ID_tri = zeros(n_np, 1);
    counter_tri = 0;
    for ny = 2 : n_np_y - 1
        for nx = 2 : n_np_x - 1
            index = (ny-1)*n_np_x + nx;
            counter_tri = counter_tri + 1;
            ID_tri(index) = counter_tri;
        end
    end
    n_eq_tri = counter_tri;% 三角形单元的自由度数
    % Triangle elements（三角形代码）
    K_tri = spalloc(n_eq_tri, n_eq_tri, 9 * n_eq_tri); % 初始化稀疏刚度矩阵
    F_tri = zeros(n_eq_tri, 1);% 初始化载荷向量

    for ee = 1 : 2*n_el
        x_ele = x_coor( IEN_tri(ee, 1:3) );
        y_ele = y_coor( IEN_tri(ee, 1:3) );
        k_ele = zeros(3, 3);% 初始化当前元素的刚度矩阵
        f_ele = zeros(3, 1);% 初始化当前元素的载荷向量
        for ll = 1 : n_int
            x_l = 0.0; y_l = 0.0;
            dx_dxi = 0.0; dx_deta = 0.0;
            dy_dxi = 0.0; dy_deta = 0.0;
            for aa = 1 : 3
                x_l = x_l + x_ele(aa) * Quad_tri(aa, xi(ll), eta(ll));
                y_l = y_l + y_ele(aa) * Quad_tri(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad_tri(aa, xi(ll), eta(ll));
                dx_dxi = dx_dxi + x_ele(aa) * Na_xi;
                dx_deta = dx_deta + x_ele(aa) * Na_eta;
                dy_dxi = dy_dxi + y_ele(aa) * Na_xi;
                dy_deta = dy_deta + y_ele(aa) * Na_eta;
            end
            detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
            for aa = 1 : 3
                Na = Quad_tri(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad_tri(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
                for bb = 1 : 3
                    Nb = Quad_tri(bb, xi(ll), eta(ll));
                    [Nb_xi, Nb_eta] = Quad_grad_tri(bb, xi(ll), eta(ll));
                    Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
                    k_ele(aa, bb) = k_ele(aa, bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
                end
            end
        end
        for aa = 1 : 3
            PP = ID_tri(IEN_tri(ee, aa));
            if PP > 0
                F_tri(PP) = F_tri(PP) + f_ele(aa);
                for bb = 1 : 3
                    QQ = ID_tri(IEN_tri(ee, bb));
                    if QQ > 0
                        K_tri(PP, QQ) = K_tri(PP, QQ) + k_ele(aa, bb);
                    end
                end
            end
        end
    end

    % Solve the stiffness matrix for triangle elements(三角形代码)
    dn_tri = K_tri \ F_tri;
    disp_tri = zeros(n_np, 1);
    for ii = 1 : n_np
        index = ID_tri(ii);
        if index > 0
            disp_tri(ii) = dn_tri(index);
        end
    end
    % % 误差函数（不用注释）
    error_function_L2 = @(u_num, u_man) sqrt(sum((u_num - u_man).^2));
    error_function_H1 = @(u_num, u_man, u_num_x, u_man_x, u_num_y, u_man_y) ...
    sqrt(sum((u_num - u_man).^2) + sum((u_num_x - u_man_x).^2) + sum((u_num_y - u_man_y).^2));
    % Error for triangle elements（三角形代码）
    u_man = manufactured_solution(x_coor, y_coor);
    u_man_x = manufact_solution_grad_x(x_coor, y_coor);
    u_man_y = manufact_solution_grad_y(x_coor, y_coor);

    errors_L2_tri(i) = error_function_L2(disp_tri, u_man);
    errors_H1_tri(i) = error_function_H1(disp_tri, u_man, ...
        exact_x(x_coor, y_coor), exact_y(x_coor, y_coor), ...
        u_man_x, u_man_y);
end

figure;
mesh_sizes = 1 ./ n_el_x_values;
% %三角形图
loglog(mesh_sizes, errors_L2_tri, '-o', 'DisplayName', 'L2 Error (Triangle)');
hold on;
loglog(mesh_sizes, errors_H1_tri, '-x', 'DisplayName', 'H1 Error (Triangle)');
slopeL2=errors_L2_tri/mesh_sizes;
slopeH1=errors_H1_tri/mesh_sizes;

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
