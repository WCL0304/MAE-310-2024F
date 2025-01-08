
% 计算解析解
function [u, v] = exact_solution_helper(x, y, L, R, type)
    % 计算解析解
    T = 10e3; % 牵引力
    r = sqrt(x.^2 + y.^2);
    theta = atan2(y, x);
    if strcmp(type, 'u')
        u = T * (1 - R^2 / r^2) * cos(theta);
        v = 0; % 只计算 u 方向的位移
    elseif strcmp(type, 'v')
        u = 0; % 只计算 v 方向的位移
        v = T * (1 - R^2 / r^2) * sin(theta);
    end
end