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

