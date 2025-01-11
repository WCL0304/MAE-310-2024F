% 自由度编号
function [n_np, ID, n_eq] = assign_dof(x_coor)
    n_np = size(x_coor, 1); % 总节点数
    ID = zeros(n_np, 2); % 每个节点两个自由度 (u, v)
    counter = 0;
    for i = 1:n_np
        counter = counter + 1;
        ID(i, 1) = counter;   % u 自由度
        counter = counter + 1;
        ID(i, 2) = counter;   % v 自由度
    end
    n_eq = max(ID(:)); % 总方程数
end