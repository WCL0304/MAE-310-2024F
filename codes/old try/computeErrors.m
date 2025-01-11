function [errors_L2, errors_H1] = computeErrors(nodes, elements, u_full, manufactured_solution, manufact_solution_grad_x, manufact_solution_grad_y)
    n_nodes = size(nodes, 1);
    n_elements = size(elements, 1);
    
    % 计算制造解的位移和梯度
    u_man = manufactured_solution(nodes(:, 1), nodes(:, 2));
    u_man_x = manufact_solution_grad_x(nodes(:, 1), nodes(:, 2));
    u_man_y = manufact_solution_grad_y(nodes(:, 1), nodes(:, 2));
    
    % 计算 L2 和 H1 误差
    errors_L2 = error_function_L2(u_full, u_man);
    errors_H1 = error_function_H1(u_full, u_man, u_man_x, u_man_x, u_man_y, u_man_y);
    
    function errors = error_function_L2(u_num, u_man)
        errors = sqrt(sum((u_num - u_man).^2));
    end
    
    function errors = error_function_H1(u_num, u_man, u_num_x, u_man_x, u_num_y, u_man_y)
        errors = sqrt(sum((u_num - u_man).^2) + sum((u_num_x - u_man_x).^2) + sum((u_num_y - u_man_y).^2));
    end
end