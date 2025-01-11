% 处理边界条件
function [dirichletBC, neumannBC] = apply_boundary_conditions(boundaryEdges, nodes, L)
    dirichletBC = [];
    neumannBC = [];
    nNodes = size(nodes, 1);
    disp(['节点数组大小:', num2str(nNodes)]); % 调试信息
    disp(['边界元素数组大小:', num2str(size(boundaryEdges, 1)), 'x', num2str(size(boundaryEdges, 2))]); % 调试信息
    disp('边界元素内容:');
    disp(boundaryEdges); % 打印边界元素的具体内容
    for i = 1:size(boundaryEdges, 1)
        edge = boundaryEdges(i, :);
        disp(['处理边界元素:', num2str(edge(1)), ' ', num2str(edge(2))]); % 调试信息
        if edge(1) > nNodes || edge(2) > nNodes
            error(['边界节点索引超出节点数组范围: edge(1) = ', num2str(edge(1)), ', edge(2) = ', num2str(edge(2)), ', nNodes = ', num2str(nNodes)]);
        end
        x1 = nodes(edge(1), 1);
        y1 = nodes(edge(1), 2);
        x2 = nodes(edge(2), 1);
        y2 = nodes(edge(2), 2);
        % 判断边界类型
        if x1 == -L / 2 || x2 == -L / 2 % 左边界
            dirichletBC = [dirichletBC; edge(1), 1, 0; edge(2), 1, 0]; % u = 0
        elseif x1 == L / 2 || x2 == L / 2 % 右边界
            neumannBC = [neumannBC; edge(1), 1, 10e3; edge(2), 1, 10e3]; % 牵引力
        elseif y1 == 0 || y2 == 0 % 底边界
            dirichletBC = [dirichletBC; edge(1), 2, 0; edge(2), 2, 0]; % v = 0
        elseif y1 == L / 2 || y2 == L / 2 % 顶边界
            % 可以根据需要添加其他边界条件
        end
    end
end