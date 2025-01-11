function [L2_error, H1_error] = calculate_errors(nodes, displacements, exact_solution, element_connectivity)
    % exact_solution 应该是一个接受 x, y 坐标和一个标识符（'u' 或 'v'）的函数，用来返回精确的位移
    % element_connectivity 是元素连通性矩阵，用于确定每个元素的节点

    % 初始化误差向量
    u_exact = zeros(size(nodes, 1), 1);
    v_exact = zeros(size(nodes, 1), 1);
    u_computed = displacements(1:2:end);
    v_computed = displacements(2:2:end);

    % 计算每个节点上的精确位移
    for i = 1:size(nodes, 1)
        x = nodes(i, 1);
        y = nodes(i, 2);
        [u_exact(i), ~] = exact_solution(x, y, 'u');
        [~, v_exact(i)] = exact_solution(x, y, 'v');
    end

    % 计算 L2 范数误差
    L2_error = sqrt(sum((u_computed - u_exact) .^ 2 + (v_computed - v_exact) .^ 2)) / sqrt(sum(u_exact .^ 2 + v_exact .^ 2));

    % 初始化 H1 范数误差的累加器
    H1_error_squared = 0;
    exact_H1_squared = 0;

    % 遍历每个元素，计算其内部的梯度误差
    for ei = 1:size(element_connectivity, 1)
        % 获取元素节点
        node_ids = element_connectivity(ei, :);
        element_nodes = nodes(node_ids, :);
        element_displacements = displacements(2*node_ids-1);
        element_displacements = [element_displacements, displacements(2*node_ids)];

        % 使用数值微分或有限元方法计算梯度
        % 这里假设用户会实现一个函数来计算梯度
        [grad_u_computed, grad_v_computed] = compute_gradient(element_nodes, element_displacements);
        
        % 计算每个节点的精确梯度
        grad_u_exact = zeros(size(element_nodes, 1), 2);
        grad_v_exact = zeros(size(element_nodes, 1), 2);
        for i = 1:size(element_nodes, 1)
            x = element_nodes(i, 1);
            y = element_nodes(i, 2);
            [grad_u_exact(i,:), ~] = exact_solution(x, y, 'grad_u');
            [~, grad_v_exact(i,:)] = exact_solution(x, y, 'grad_v');
        end

        % 在每个元素内部积分梯度误差
        for i = 1:size(element_nodes, 1)
            H1_error_squared = H1_error_squared + (norm(grad_u_computed(i,:) - grad_u_exact(i,:)))^2 + (norm(grad_v_computed(i,:) - grad_v_exact(i,:)))^2;
            exact_H1_squared = exact_H1_squared + (norm(grad_u_exact(i,:)))^2 + (norm(grad_v_exact(i,:)))^2;
        end
    end

    % 计算最终的 H1 范数误差
    H1_error = sqrt(H1_error_squared) / sqrt(exact_H1_squared + eps); % 添加 eps 防止分母为零
end

% 假设的函数：计算梯度，用户需要实现此部分
function [grad_u, grad_v] = compute_gradient(nodes, displacements)
    % 实现计算节点梯度的逻辑
    % 可以使用有限差分、有限元方法中的局部梯度计算等
    % 这里只是一个示例
    % 实际应用中可能需要根据具体的单元类型（如线性三角形或双线性四边形）调整
    grad_u = zeros(size(nodes, 1), 2);
    grad_v = zeros(size(nodes, 1), 2);
end