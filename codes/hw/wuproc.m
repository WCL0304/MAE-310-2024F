clc;
clear;
close all;
%不知道为什么画不出直线

elements = [2, 4, 6, 8, 10, 12, 14, 16]; 
orders = [1, 2, 3];% 定义不同的元素阶次：一次、二次、三次
h_values = 1 ./ elements; 

for p = orders
    errors_L2 = zeros(size(elements));% 初始化L2误差数组
    errors_H1 = zeros(size(elements));% 初始化H1误差数组
    fprintf('元素阶次: %d\n', p);% 打印当前元素阶次
    for i = 1:length(elements)
        n = elements(i);
        [u_h, u_h_x, u_exact, u_exact_x] = solveFEM(n, p); % 求解
        errors_L2(i) = computeErrorL2(u_h, u_exact);% 计算L2误差
        errors_H1(i) = computeErrorH1(u_h_x, u_exact_x);% 计算H1误差
    end

    if any(errors_L2 <= 0) || any(errors_H1 <= 0)% 检查误差是否为正数
        error('Error values must be positive and non-zero.');
    end
    
    % 绘制 log(error)-log(h) 图
    figure;
    loglog(h_values, errors_L2, '-o', 'DisplayName', 'L2 Error');
    hold on;
    loglog(h_values, errors_H1, '-x', 'DisplayName', 'H1 Error');
    xlabel('log(h)');
    ylabel('log(Error)');
    legend;
    title(['误差分析 (元素阶次: ', num2str(p), ')']);
    grid on;
    
    % 计算斜率
    slope_L2 = polyfit(log(h_values), log(errors_L2), 1);
    slope_H1 = polyfit(log(h_values), log(errors_H1), 1);
    
    fprintf('L2 误差斜率: %.2f\n', slope_L2(1));
    fprintf('H1 误差斜率: %.2f\n', slope_H1(1));
end

% 定义solveFEM函数
function [u_h, u_h_x, u_exact, u_exact_x] = solveFEM(n, p)
    if p == 1
        x_nodes = linspace(0, 1, n + 1);
    else
        % 对于高阶元素，添加中间节点
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
    u_nodes = x_nodes.^2; % 精确解值
    
    % 数值解插值函数
    u_h = @(x) interp1(x_nodes, u_nodes, x, interpolationMethod(p), 'extrap'); % 数值解插值
    
    % 计算数值解导数
    u_h_x = @(x) numericalDerivative(u_h, x);
    
    u_exact = @(x) x.^2; % 精确解
    u_exact_x = @(x) 2 * x; % 精确导数
end

% 定义interpolationMethod函数
function method = interpolationMethod(p)
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

% 定义numericalDerivative函数
function value = numericalDerivative(u_h, x)
    % 使用五点差分计算导数
    h = 1e-4;
    value = (-u_h(x + 2*h) + 8*u_h(x + h) - 8*u_h(x - h) + u_h(x - 2*h)) / (12*h);
end

% 定义computeErrorL2函数
function eL2 = computeErrorL2(u_h, u_exact)
    numerator = Gauss(@(x) (u_h(x) - u_exact(x)).^2);
    denominator = Gauss(@(x) (u_exact(x)).^2);

    if denominator == 0
        error('Denominator in L2 error calculation is zero.');
    end

    eL2 = sqrt(numerator / denominator);
end

% 定义computeErrorH1函数
function eH1 = computeErrorH1(u_h_x, u_exact_x)
    numerator = Gauss(@(x) (u_h_x(x) - u_exact_x(x)).^2);
    denominator = Gauss(@(x) (u_exact_x(x)).^2);

    if denominator == 0
        error('Denominator in H1 error calculation is zero.');
    end

    eH1 = sqrt(numerator / denominator);
end

% Gauss.m 实现
function integral = Gauss(func)% 高斯积分
    % 对 [0, 1] 区间进行积分
    [points, weights] = gaussPointsAndWeights(10); % 使用 10 点高斯积分
    points = (points + 1) / 2; % 映射到 [0, 1]
    weights = weights / 2;    % 权重调整
    integral = sum(weights .* func(points));
end

% gaussPointsAndWeights.m 实现
function [points, weights] = gaussPointsAndWeights(n)
    if n == 10
        points = [-0.9739065285, -0.8650633666, -0.6794095682, -0.4333953941, -0.1488743390, ...
                   0.1488743390,  0.4333953941,  0.6794095682,  0.8650633666,  0.9739065285];
        weights = [0.0666713443,  0.1494513491,  0.2190863625,  0.2692667193,  0.2955242247, ...
                   0.2955242247,  0.2692667193,  0.2190863625,  0.1494513491,  0.0666713443];
    else
        error('高斯积分点数未实现');
    end
end
