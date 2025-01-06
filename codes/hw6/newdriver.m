clear;
clc;
close all;

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation
n_en   = 4;               % number of nodes in an element
n_el_x = 60;               % number of elements in x-dir
n_el_y = 60;               % number of elements in y-dir
n_el   = n_el_x * n_el_y; % total number of elements

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% generate the nodal coordinates
for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index
    x_coor(index) = (nx-1) * hx;
    y_coor(index) = (ny-1) * hy;
  end
end
% IEN array for triangles
IEN = zeros(2*n_el, 3); % 每个四边形元素被分成两个三角形元素
for ex = 1 : n_el_x
    for ey = 1 : n_el_y
        ee = (ey-1) * n_el_x + ex; % 四边形元素索引
        IEN(2*(ee-1)+1, :) = [(ey-1) * n_np_x + ex, (ey-1) * n_np_x + ex + 1, ey * n_np_x + ex + 1]; % 第一个三角形元素
        IEN(2*ee, :) = [(ey-1) * n_np_x + ex, ey * n_np_x + ex + 1, ey * n_np_x + ex]; % 第二个三角形元素
    end
end

% ID array for triangles
ID = zeros(n_np, 1);
counter = 0;
for ny = 2 : n_np_y - 1
    for nx = 2 : n_np_x - 1
        index = (ny-1)*n_np_x + nx;
        counter = counter + 1;
        ID(index) = counter;
    end
end
n_eq = counter;

LM = ID(IEN);

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : 2*n_el % 三角形元素数量
    x_ele = x_coor( IEN(ee, 1:3) );
    y_ele = y_coor( IEN(ee, 1:3) );
    k_ele = zeros(3, 3); % element stiffness matrix for triangles
    f_ele = zeros(3, 1); % element load vector for triangles
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
            end % end of bb loop
        end % end of aa loop
    end % end of quadrature loop
    for aa = 1 : 3
        PP = LM(ee, aa);
        if PP > 0
            F(PP) = F(PP) + f_ele(aa);
            for bb = 1 : 3
                QQ = LM(ee, bb);
                if QQ > 0
                    K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
                else
                    % modify F with the boundary data
                    % here we do nothing because the boundary data g is zero or
                    % homogeneous
                end
            end
        end
    end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 1);
for ii = 1 : n_np
    index = ID(ii);
    if index > 0
        disp(ii) = dn(index);
    else
        % modify disp with the g data. Here it does nothing because g is zero
    end
end

% 制造解
manufactured_solution = @(x, y) x.^2 .* y.^2;

% 误差函数
error_function = @(u_num, u_man) sqrt(sum(sum((u_num - u_man).^2)));

% 初始网格尺寸
n_el_x_values = [10, 20, 40, 80, 160];
n_el_y_values = n_el_x_values;

% 存储误差
errors_L2_quad = zeros(size(n_el_x_values));
errors_L2_tri = zeros(size(n_el_x_values));

% 计算不同网格尺寸下的误差
for i = 1:length(n_el_x_values)
    n_el_x = n_el_x_values(i);
    n_el_y = n_el_y_values(i);
    n_el = n_el_x * n_el_y;
    
    n_np_x = n_el_x + 1;
    n_np_y = n_el_y + 1;
    n_np = n_np_x * n_np_y;
    
    x_coor = zeros(n_np, 1);
    y_coor = x_coor;
    hx = 1.0 / n_el_x;
    hy = 1.0 / n_el_y;
    
    % 生成节点坐标
    for ny = 1 : n_np_y
        for nx = 1 : n_np_x
            index = (ny-1)*n_np_x + nx;
            x_coor(index) = (nx-1) * hx;
            y_coor(index) = (ny-1) * hy;
        end
    end
    
    % 重新定义IEN数组和ID数组
    IEN_quad = zeros(n_el, 4);
    for ex = 1 : n_el_x
        for ey = 1 : n_el_y
            ee = (ey-1) * n_el_x + ex;
            IEN_quad(ee, :) = [(ey-1) * n_np_x + ex, (ey-1) * n_np_x + ex + 1, ey * n_np_x + ex + 1, ey * n_np_x + ex];
        end
    end
    
    IEN_tri = zeros(2*n_el, 3);
    for ex = 1 : n_el_x
        for ey = 1 : n_el_y
            ee = (ey-1) * n_el_x + ex;
            IEN_tri(2*(ee-1)+1, :) = [(ey-1) * n_np_x + ex, (ey-1) * n_np_x + ex + 1, ey * n_np_x + ex + 1];
            IEN_tri(2*ee, :) = [(ey-1) * n_np_x + ex, ey * n_np_x + ex + 1, ey * n_np_x + ex];
        end
    end
    
    ID_quad = zeros(n_np, 1);
    ID_tri = zeros(n_np, 1);
    counter_quad = 0;
    counter_tri = 0;
    for ny = 2 : n_np_y - 1
        for nx = 2 : n_np_x - 1
            index = (ny-1)*n_np_x + nx;
            counter_quad = counter_quad + 1;
            counter_tri = counter_tri + 1;
            ID_quad(index) = counter_quad;
            ID_tri(index) = counter_tri;
        end
    end
    
    n_eq_quad = counter_quad;
    n_eq_tri = counter_tri;
    
    % 组装K和F的代码
    [xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);
    
    % Quadrilateral elements
    K_quad = spalloc(n_eq_quad, n_eq_quad, 9 * n_eq_quad);
    F_quad = zeros(n_eq_quad, 1);
    
    for ee = 1 : n_el
        x_ele = x_coor( IEN_quad(ee, 1:4) );
        y_ele = y_coor( IEN_quad(ee, 1:4) );
        k_ele = zeros(4, 4);
        f_ele = zeros(4, 1);
        for ll = 1 : n_int
            x_l = 0.0; y_l = 0.0;
            dx_dxi = 0.0; dx_deta = 0.0;
            dy_dxi = 0.0; dy_deta = 0.0;
            for aa = 1 : 4
                x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
                y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                dx_dxi = dx_dxi + x_ele(aa) * Na_xi;
                dx_deta = dx_deta + x_ele(aa) * Na_eta;
                dy_dxi = dy_dxi + y_ele(aa) * Na_xi;
                dy_deta = dy_deta + y_ele(aa) * Na_eta;
            end
            detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
            for aa = 1 : 4
                Na = Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
                for bb = 1 : 4
                    Nb = Quad(bb, xi(ll), eta(ll));
                    [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                    Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
                    k_ele(aa, bb) = k_ele(aa, bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
                end
            end
        end
        for aa = 1 : 4
            PP = ID_quad(IEN_quad(ee, aa));
            if PP > 0
                F_quad(PP) = F_quad(PP) + f_ele(aa);
                for bb = 1 : 4
                    QQ = ID_quad(IEN_quad(ee, bb));
                    if QQ > 0
                        K_quad(PP, QQ) = K_quad(PP, QQ) + k_ele(aa, bb);
                    end
                end
            end
        end
    end
    
    % Solve the stiffness matrix for quadrilateral elements
    dn_quad = K_quad \ F_quad;
    disp_quad = zeros(n_np, 1);
    for ii = 1 : n_np
        index = ID_quad(ii);
        if index > 0
            disp_quad(ii) = dn_quad(index);
        end
    end
    
    % Error for quadrilateral elements
    u_man = manufactured_solution(x_coor, y_coor);
    e_quad = disp_quad - u_man;
    errors_L2_quad(i) = error_function(e_quad, 0);
    
    % Triangle elements
    K_tri = spalloc(n_eq_tri, n_eq_tri, 9 * n_eq_tri);
    F_tri = zeros(n_eq_tri, 1);
    
    for ee = 1 : 2*n_el
        x_ele = x_coor( IEN_tri(ee, 1:3) );
        y_ele = y_coor( IEN_tri(ee, 1:3) );
        k_ele = zeros(3, 3);
        f_ele = zeros(3, 1);
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
    
    % Error for triangle elements
    u_man = manufactured_solution(x_coor, y_coor);
    e_tri = disp_tri - u_man;
    errors_L2_tri(i) = error_function(e_tri, 0);
end

% 绘制误差与网格尺寸的对数-对数图
mesh_sizes = 1.0 ./ n_el_x_values;
loglog(mesh_sizes, errors_L2_quad, '-o', 'DisplayName', 'Quadrilateral Elements');
hold on;
loglog(mesh_sizes, errors_L2_tri, '-x', 'DisplayName', 'Triangle Elements');
xlabel('网格尺寸');
ylabel('L2误差');
title('误差与网格尺寸的关系');
legend;
grid on;