clc; 
clear;
close all;
%警告不知道为什么越debug越多
%结果很奇怪
elements = [10, 20, 40, 80];
orders = 1; 
h_values = 1 ./ elements; 
tol_values = [1e-2, 1e-4, 1e-6]; 

for tol = tol_values
    errors_gmres = zeros(size(elements));% 初始化GMRES误差数组
    errors_direct = zeros(size(elements));% 初始化直接法误差数组
    iterations = zeros(size(elements));% 初始化迭代次数数组
    fprintf('GMRES 方法，容差 tol = %e\n', tol);% 打印当前容差值
    for i = 1:length(elements)
        n = elements(i);
        [A, b, u_exact_func] = assembleFEM(n, orders); % 调用assembleFEM函数组装有限元问题
        u_exact = u_exact_func(linspace(0, 1, n + 1))';
        
        % GMRES 方法求解
        restart = [];
        maxit = 10000;
        [u_gmres, flag, relres, iter] = gmres(A, b, restart, tol, maxit);
        iterations(i) = sum(iter);
        
        % 直接方法求解
        u_direct = A \ b;
        
        % 计算误差
        errors_gmres(i) = norm(u_gmres - u_exact, 2) / norm(u_exact, 2);% 计算GMRES的相对L2误差
        errors_direct(i) = norm(u_direct - u_exact, 2) / norm(u_exact, 2);% 计算直接法的相对L2误差
    end

    

    for i = 1:length(elements)
        fprintf('网格数: %d, GMRES 误差: %.2e, 直接法误差: %.2e, 迭代次数: %d\n', ...
            elements(i), errors_gmres(i), errors_direct(i), iterations(i));
    end
    fprintf('\n');
end

% 定义assembleFEM函数
function [A, b, u_exact_func] = assembleFEM(n, p)

    h = 1 / n;
    x = linspace(0, 1, n + 1)';

    A = sparse(n + 1, n + 1);% 初始化系数矩阵A
    b = zeros(n + 1, 1);% 初始化向量b

    u_exact_func = @(x) x.^2 - x;
    for i = 2:n
        A(i, i - 1) = -1 / h;
        A(i, i)     = 2 / h;
        A(i, i + 1) = -1 / h;
        b(i) = 2 * h; 
    end

    A(1, 1) = 1;
    b(1) = 0;
    A(n + 1, n + 1) = 1;
    b(n + 1) = 0;
end
