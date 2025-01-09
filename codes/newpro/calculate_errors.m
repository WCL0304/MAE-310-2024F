% 计算误差
function [L2_error, H1_error] = calculate_errors(nodes, displacements, exact_solution)
    u = displacements(1:2:end);
    v = displacements(2:2:end);
    
    % 计算解析解
    [u_exact, v_exact] = deal(arrayfun(@(x, y) exact_solution(x, y, 'u'), nodes(:, 1), nodes(:, 2)));
    [u_0, v_0] = deal(arrayfun(@(x, y) exact_solution(x, y, 'v'), nodes(:, 1), nodes(:, 2)));
    
    % L2范数误差
    L2_error_u = norm(u - u_exact) / norm(u_exact);
    L2_error_v = norm(v - v_0) / norm(v_0);
    L2_error = mean([L2_error_u, L2_error_v]);
    
    % 计算梯度
    [dudx, dudy] = gradient(u, nodes(:, 1), nodes(:, 2));
    [dvdx, dvdy] = gradient(v, nodes(:, 1), nodes(:, 2));
    [dudx_exact, dudy_exact] = gradient(u_exact, nodes(:, 1), nodes(:, 2));
    [dvdx_exact, dvdy_exact] = gradient(v_0, nodes(:, 1), nodes(:, 2));
    
    % H1范数误差
    H1_error_u = norm([dudx - dudx_exact; dudy - dudy_exact]) / norm([dudx_exact; dudy_exact]);
    H1_error_v = norm([dvdx - dvdx_exact; dvdy - dvdy_exact]) / norm([dvdx_exact; dvdy_exact]);
    H1_error = mean([H1_error_u, H1_error_v]);
end