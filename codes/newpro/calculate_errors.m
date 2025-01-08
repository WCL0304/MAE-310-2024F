function [L2_error, H1_error] = calculate_errors(nodes, displacements, exact_solution)
    % 计算 L2 和 H1 范数误差
    nNodes = size(nodes, 1);
    u = displacements(:, 1);
    v = displacements(:, 2);
    u_exact = exact_solution(nodes, 'u');
    v_exact = exact_solution(nodes, 'v');

    % L2 范数误差
    L2_error_u = norm(u - u_exact) / norm(u_exact);
    L2_error_v = norm(v - v_exact) / norm(v_exact);
    L2_error = mean([L2_error_u, L2_error_v]);

    % H1 范数误差
    grad_u = gradient(u, nodes(:, 1), nodes(:, 2));
    grad_v = gradient(v, nodes(:, 1), nodes(:, 2));
    grad_u_exact = gradient(u_exact, nodes(:, 1), nodes(:, 2));
    grad_v_exact = gradient(v_exact, nodes(:, 1), nodes(:, 2));

    H1_error_u = norm([grad_u(1) - grad_u_exact(1); grad_u(2) - grad_u_exact(2)]) / norm(grad_u_exact);
    H1_error_v = norm([grad_v(1) - grad_v_exact(1); grad_v(2) - grad_v_exact(2)]) / norm(grad_v_exact);
    H1_error = mean([H1_error_u, H1_error_v]);
end