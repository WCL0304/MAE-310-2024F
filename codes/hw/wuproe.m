clc; 
clear;
close all;
%不知道是高斯积分点算错了，还是误差算错了，出的图很奇怪
elements = [2, 4, 6, 8, 10, 12, 14, 16]; % 网格数
orders = [1, 2, 3]; % 线性、二次、三次
h_values = 1 ./ elements; % 网格大小

gauss_points = [1, 2, 3, 4, 5, 6]; % 高斯积分点数
fprintf('实验不同的高斯积分点数：\n');
for n_gauss = gauss_points
    fprintf('积分点数: %d\n', n_gauss);
    errors_L2 = zeros(size(elements));
    for i = 1:length(elements)
        n = elements(i);
        [u_h, u_h_x, u_exact, u_exact_x] = solveFEM(n, 3); % 三次元素
        errors_L2(i) = computeErrorL2(@(x) u_h(x), @(x) u_exact(x), n_gauss);
    end
    
    figure;
    loglog(h_values, errors_L2, '-o');
    title(['高斯积分点数: ', num2str(n_gauss)]);
    xlabel('h');
    ylabel('L2 Error');
    grid on;
end

fprintf('GMRES与直接求解对比：\n');
tol = 1e-6;
for i = 1:length(elements)
    n = elements(i);
    [A, b] = constructSystem(n, 3); % 构造有限元系统
    max_iter = size(A, 1); % 最大迭代次数
    % 直接解法
    u_direct = A \ b;
    % GMRES解法
    [u_gmres, flag] = gmres(A, b, [], tol, max_iter);
    fprintf('网格数: %d, GMRES flag: %d\n', n, flag);
end


% solveFEM: 构造有限元解
function [u_h, u_h_x, u_exact, u_exact_x] = solveFEM(n, p)
    if p == 1
        x_nodes = linspace(0, 1, n + 1);
    else
        % 高阶元素，添加中间节点
        x_elements = linspace(0, 1, n + 1);
        xi = linspace(0, 1, p + 1);
        x_nodes = [];
        for i = 1:n
            x_e = x_elements(i) + (x_elements(i+1) - x_elements(i)) * xi;
            x_nodes = [x_nodes, x_e(1:end-1)];
        end
        x_nodes = [x_nodes, x_elements(end)];
        x_nodes = unique(x_nodes); % 去除重复节点
    end
    u_nodes = x_nodes.^2; % 节点处的精确解值
    
    % 数值解插值函数
    u_h = @(x) interp1(x_nodes, u_nodes, x, interpolationMethod(p), 'extrap'); % 数值解插值
    
    % 计算数值解的导数
    u_h_x = @(x) numericalDerivative(u_h, x);
    
    u_exact = @(x) x.^2; % 精确解
    u_exact_x = @(x) 2 * x; % 精确导数
end

% interpolationMethod: 插值
function method = interpolationMethod(p)
    % 返回插值方法
    if p == 1
        method = 'linear';
    elseif p == 2
        method = 'pchip'; % 分段三次 Hermite 插值
    elseif p == 3
        method = 'spline'; % 三次样条插值
    else
        error('Unsupported element order');
    end
end

% numericalDerivative: 数值微分
function value = numericalDerivative(u_h, x)
    % 使用五点差分计算导数
    h = 1e-4;
    value = (-u_h(x + 2*h) + 8*u_h(x + h) - 8*u_h(x - h) + u_h(x - 2*h)) / (12*h);
end

function eL2 = computeErrorL2(u_h, u_exact, n_gauss)
    numerator = Gauss(@(x) (u_h(x) - u_exact(x)).^2, n_gauss);
    denominator = Gauss(@(x) (u_exact(x)).^2, n_gauss);

    if denominator == 0
        error('Denominator in L2 error calculation is zero.');
    end

    eL2 = sqrt(numerator / denominator);
end

% constructSystem: 构造有限元
function [A, b] = constructSystem(n, p)
    A = diag(2*ones(1, n+1)) - diag(ones(1, n), 1) - diag(ones(1, n), -1); 
    b = ones(n+1, 1); 
end

% Gauss: 高斯积分实现
function integral = Gauss(func, n)
    [points, weights] = gaussPointsAndWeights(n); % 获取积分点和权重
    points = (points + 1) / 2; % 映射到 [0, 1]
    weights = weights / 2;    % 权重调整
    integral = sum(weights .* func(points));
end

% gaussPointsAndWeights: 返回高斯积分点和权重
function [points, weights] = gaussPointsAndWeights(n)
    if n == 1
        points = 0;
        weights = 2;
    elseif n == 2
        points = [-0.5773502692, 0.5773502692];
        weights = [1, 1];
    elseif n == 3
        points = [-0.7745966692, 0, 0.7745966692];
        weights = [0.5555555556, 0.8888888889, 0.5555555556];
    elseif n == 4
        points = [-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116];
        weights = [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451];
    elseif n == 5
        points = [-0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459];
        weights = [0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851];
    elseif n == 6
        points = [-0.9324695142, -0.6612093865, -0.2386191861, 0.2386191861, 0.6612093865, 0.9324695142];
        weights = [0.1713244924, 0.3607615730, 0.4679139346, 0.4679139346, 0.3607615730, 0.1713244924];
    else
        error('Unsupported number of Gauss points');
    end
end
