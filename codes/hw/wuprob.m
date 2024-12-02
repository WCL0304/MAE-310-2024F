clc;
clear;
close all;
%driver中完成的大部分内容不重复，这里就画图
% 问题定义
elements = [2, 4, 6, 8, 10, 12, 14, 16];% 定义一系列元素数量，用于不同网格尺寸的测试
errors_L2 = zeros(size(elements));% 初始化一个数组，用于存储每个网格尺寸下的L2误差
errors_H1 = zeros(size(elements));% 初始化一个数组，用于存储每个网格尺寸下的H1误差
h_values = 1 ./ elements;% 计算每个元素数量对应的网格尺寸h

for i = 1:length(elements)% 遍历每个元素数量
    n = elements(i);% 获取当前元素数量
    [u_h, u_h_x, u_exact, u_exact_x, x_nodes] = solveFEM(n);% 调用solveFEM函数求解有限元问题
    errors_L2(i) = computeErrorL2(u_h, u_exact, x_nodes);% 计算并存储L2误差
    errors_H1(i) = computeErrorH1(u_h_x, u_exact_x, x_nodes);% 计算并存储H1误差
end

% 绘制 log(error)-log(h) 图
figure;

% 绘制 L2 误差
subplot(2,1,1);
loglog(h_values, errors_L2, '-o', 'DisplayName', 'L2 Error');% 绘制L2误差的对数图
xlabel('log(h)');
ylabel('log(L2 Error)');
title('L2 Error vs. Mesh Size');
legend show;

% 绘制 H1 误差
subplot(2,1,2);
loglog(h_values, errors_H1, '-x', 'DisplayName', 'H1 Error');% 绘制H1误差的对数图
xlabel('log(h)');
ylabel('log(H1 Error)');
title('H1 Error vs. Mesh Size');
legend show;


slope_L2 = polyfit(log(h_values), log(errors_L2), 1);
slope_H1 = polyfit(log(h_values), log(errors_H1), 1);

fprintf('L2 Error Slope: %.2f\n', slope_L2(1));
fprintf('H1 Error Slope: %.2f\n', slope_H1(1));

% 定义solveFEM函数
function [u_h, u_h_x, u_exact, u_exact_x, x_nodes] = solveFEM(n)

    x_nodes = linspace(0, 1, n + 1); % 在区间[0,1]上均匀划分n+1个节点
    u_nodes = x_nodes.^2; % 计算每个节点处的精确解的值
    u_h = @(x) interp1(x_nodes, u_nodes, x, 'linear'); % 构造一个线性插值函数，作为数值解

    slopes = diff(u_nodes) ./ diff(x_nodes); % 计算每个元素上的斜率
    u_h_x = @(x) piecewiseDerivative(x, x_nodes, slopes);% 构造一个分段导数函数，作为数值解的导数

    u_exact = @(x) x.^2; % 定义精确解函数
    u_exact_x = @(x) 2 * x; % 定义精确解的导数函数
end
% 定义piecewiseDerivative函数：计算分段导数
function value = piecewiseDerivative(x, x_nodes, slopes)% 初始化输出数组，其大小与输入数组x相同
    
    value = zeros(size(x));
    for i = 1:length(x)
        xi = x(i);
        idx = find(xi >= x_nodes(1:end-1) & xi <= x_nodes(2:end));
        if ~isempty(idx)% 如果当前点在某个元素内
            value(i) = slopes(idx);
        else
            if xi < x_nodes(1) % 如果当前点小于第一个节点
                value(i) = slopes(1);
            elseif xi > x_nodes(end)% 如果当前点大于最后一个节点
                value(i) = slopes(end);
            end
        end
    end
end

% 定义computeErrorL2函数：计算L2误差
function eL2 = computeErrorL2(u_h, u_exact, x_nodes)
    e_square = 0;% 初始化误差平方和为0
    for i = 1:length(x_nodes)-1
        a = x_nodes(i);
        b = x_nodes(i+1);
        e_square = e_square + GaussLegendre5(@(x) (u_h(x) - u_exact(x)).^2, a, b);% 计算当前元素上的误差平方，并累加到总误差平方和中
    end
    eL2 = sqrt(e_square);% 计算总误差的平方根，得到L2误差
end
% 定义computeErrorH1函数：计算H1误差
function eH1 = computeErrorH1(u_h_x, u_exact_x, x_nodes)
    e_square = 0;% 初始化误差平方和为0
    for i = 1:length(x_nodes)-1
        a = x_nodes(i);
        b = x_nodes(i+1);
        e_square = e_square + GaussLegendre5(@(x) (u_h_x(x) - u_exact_x(x)).^2, a, b);% 计算当前元素上的误差平方，并累加到总误差平方和中
    end
    eH1 = sqrt(e_square);% 计算总误差的平方根，得到H1误差
end


function integral = GaussLegendre5(func, a, b)
    [xi, weights] = GaussLegendre(5);
    xi = (b - a) / 2 * xi + (b + a) / 2;
    weights = (b - a) / 2 * weights;
    integral = sum(weights .* arrayfun(func, xi));
end

function [xi, weights] = GaussLegendre(n)
    switch n
        case 5
            xi = [-0.906179845938664, -0.538469310105683, 0, ...
                   0.538469310105683, 0.906179845938664];
            weights = [0.236926885056189, 0.478628670499366, 0.568888888888889, ...
                       0.478628670499366, 0.236926885056189];
        otherwise
            error('Only 5-point Gauss-Legendre implemented.');
    end
end


function value = PolyShape(degree, node, xi, derivative)

    if degree == 1 
        if derivative == 0  
            if node == 1
                value = (1 - xi) / 2;
            elseif node == 2
                value = (1 + xi) / 2;
            else
                error('Invalid node index for linear element.');
            end
        elseif derivative == 1  
            if node == 1
                value = -0.5;
            elseif node == 2
                value = 0.5;
            else
                error('Invalid node index for linear element.');
            end
        end
    else
        error('Only linear elements (degree=1) are implemented.');
    end
end
