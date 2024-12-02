clc;
clear;
close all;
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