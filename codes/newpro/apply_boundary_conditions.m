% 处理边界条件
function [dirichletBC, neumannBC] = apply_boundary_conditions(boundaryEdges, nodes, L)
    dirichletBC = [];
    neumannBC = [];
    traction_force = 10e3; % 牵引力, N/m
    
    for i = 1:size(boundaryEdges, 1)
        [x1, y1] = nodes(boundaryEdges(i, 1), :);
        [x2, y2] = nodes(boundaryEdges(i, 2), :);
        
        % 判断边界类型
        if x1 == -L / 2 || x2 == -L / 2 % 左边界
            dirichletBC = [dirichletBC; boundaryEdges(i, 1), 1, 0; boundaryEdges(i, 2), 1, 0]; % u=0
        elseif x1 == L / 2 || x2 == L / 2 % 右边界
            neumannBC = [neumannBC; boundaryEdges(i, 1), 1, traction_force; boundaryEdges(i, 2), 1, traction_force]; % 牵引力
        elseif y1 == 0 || y2 == 0 % 底边界
            dirichletBC = [dirichletBC; boundaryEdges(i, 1), 2, 0; boundaryEdges(i, 2), 2, 0]; % v=0
        elseif y1 == L / 2 || y2 == L / 2 % 顶边界
            % 顶边界无需边界条件
        else
            error(['未知边界类型: edge(1)=%d, edge(2)=%d，坐标: (%f, %f) 和 (%f, %f)', boundaryEdges(i, 1), boundaryEdges(i, 2), x1, y1, x2, y2]);
        end
    end
end