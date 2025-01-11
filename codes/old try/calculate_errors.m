% 计算L2范数误差和H1范数误差
function [L2_error, H1_error] = calculate_errors(nodes, displacements, exact_solution)
    u_exact = zeros(size(nodes, 1), 1);
    v_exact = zeros(size(nodes, 1), 1);
    u_computed = displacements(1:2:end);
    v_computed = displacements(2:2:end);

    for i = 1:size(nodes, 1)
        x = nodes(i, 1);
        y = nodes(i, 2);
        [u_exact(i), ~] = exact_solution(x, y, 'u');
        [~, v_exact(i)] = exact_solution(x, y, 'v');
    end

    L2_error = sqrt(sum((u_computed - u_exact) .^ 2 + (v_computed - v_exact) .^ 2)) / sqrt(sum(u_exact .^ 2 + v_exact .^ 2));
    H1_error = sqrt(sum((u_computed - u_exact) .^ 2 + (v_computed - v_exact) .^ 2 + ((u_computed(2:end) - u_computed(1:end-1)) - (u_exact(2:end) - u_exact(1:end-1))) .^ 2 + ((v_computed(2:end) - v_computed(1:end-1)) - (v_exact(2:end) - v_exact(1:end-1))) .^ 2)) / sqrt(sum(u_exact .^ 2 + v_exact .^ 2 + (u_exact(2:end) - u_exact(1:end-1)) .^ 2 + (v_exact(2:end) - v_exact(1:end-1)) .^ 2);
end