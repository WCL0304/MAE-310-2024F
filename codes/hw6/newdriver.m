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

