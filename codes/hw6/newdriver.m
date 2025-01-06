clear;
clc;
close all;

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x.*(1-x).*y.*(1-y);
exact_x = @(x,y) (1-2*x).*y.*(1-y);
exact_y = @(x,y) x.*(1-x).*(1-2*y);

f = @(x,y) 2.0*kappa*x.*(1-x) + 2.0*kappa*y.*(1-y); % source term

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation
n_el_x_values = [10, 20, 40, 80, 160];
n_el_y_values = n_el_x_values;

% 存储误差
errors_L2_quad = zeros(size(n_el_x_values));
errors_H1_quad = zeros(size(n_el_x_values));
errors_L2_tri = zeros(size(n_el_x_values));
errors_H1_tri = zeros(size(n_el_x_values));

% 制造的解
manufactured_solution = @(x, y) x.^2 .* y.^2;

% 制造的解的梯度
manufact_solution_grad_x = @(x, y) 2 * x .* y.^2;
manufact_solution_grad_y = @(x, y) 2 * y .* x.^2;

% 计算不同网格尺寸下的误差
for i = 1:length(n_el_x_values)
    n_el_x = n_el_x_values(i);% 当前网格的 x 方向元素数
    n_el_y = n_el_y_values(i);% 当前网格的 y 方向元素数
    n_el = n_el_x * n_el_y;
    
    n_np_x = n_el_x + 1;% x 方向的节点数
    n_np_y = n_el_y + 1;% y 方向的节点数
    n_np = n_np_x * n_np_y;
    
    hx = 1.0 / n_el_x;
    hy = 1.0 / n_el_y;
    
    % 生成节点坐标
    x_coor = zeros(n_np, 1);
    y_coor = x_coor;
    for ny = 1 : n_np_y
        for nx = 1 : n_np_x
            index = (ny-1)*n_np_x + nx;
            x_coor(index) = (nx-1) * hx;
            y_coor(index) = (ny-1) * hy;
        end
    end
    
    % 重新定义IEN数组和ID数组
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

    % Triangle elements
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

        % Solve the stiffness matrix for triangle elements
    dn_tri = K_tri \ F_tri;
    disp_tri = zeros(n_np, 1);
    for ii = 1 : n_np
        index = ID_tri(ii);
        if index > 0
            disp_tri(ii) = dn_tri(index);
        end
    end
    % 误差函数
    error_function_L2 = @(u_num, u_man) sqrt(sum((u_num - u_man).^2));
    error_function_H1 = @(u_num, u_man, u_num_x, u_man_x, u_num_y, u_man_y) ...
    sqrt(sum((u_num - u_man).^2) + sum((u_num_x - u_man_x).^2) + sum((u_num_y - u_man_y).^2));
    % Error for triangle elements
    u_man = manufactured_solution(x_coor, y_coor);
    u_man_x = manufact_solution_grad_x(x_coor, y_coor);
    u_man_y = manufact_solution_grad_y(x_coor, y_coor);
    
    errors_L2_tri(i) = error_function_L2(disp_tri, u_man);
    errors_H1_tri(i) = error_function_H1(disp_tri, u_man, ...
        exact_x(x_coor, y_coor), exact_y(x_coor, y_coor), ...
        u_man_x, u_man_y);
end

% 绘制log-log误差图
figure;
mesh_sizes = 1 ./ n_el_x_values;
loglog(mesh_sizes, errors_L2_tri, '-o', 'DisplayName', 'L2 Error (Triangle)');
hold on;
loglog(mesh_sizes, errors_H1_tri, '-x', 'DisplayName', 'H1 Error (Triangle)');
xlabel('Mesh Size (h)');
ylabel('Error');
title('Convergence Rate of Errors');
legend;
grid on;

    