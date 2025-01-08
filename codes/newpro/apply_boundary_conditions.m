function [dirichletBC, neumannBC] = apply_boundary_conditions(boundaryEdges, nodes, boundaryConditions)
    % 处理边界条件
    dirichletBC = [];
    neumannBC = [];

    for i = 1:size(boundaryEdges, 1)
        edge = boundaryEdges(i, :);
        x1 = nodes(edge(1), 1);
        y1 = nodes(edge(1), 2);
        x2 = nodes(edge(2), 1);
        y2 = nodes(edge(2), 2);

        % 判断边界类型
        if x1 == -L/2 || x2 == -L/2 % 左边界
            dirichletBC = [dirichletBC; edge(1), 1, 0; edge(2), 1, 0]; % u = 0
        elseif x1 == L/2 || x2 == L/2 % 右边界
            neumannBC = [neumannBC; edge(1), 1, 10e3; edge(2), 1, 10e3]; % 牵引力
        elseif y1 == 0 || y2 == 0 % 底边界
            dirichletBC = [dirichletBC; edge(1), 2, 0; edge(2), 2, 0]; % v = 0
        elseif y1 == L/2 || y2 == L/2 % 顶边界
            dirichletBC = [dirichletBC; edge(1), 2, 0; edge(2), 2, 0]; % v = 0
        end
    end
end