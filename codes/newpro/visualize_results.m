% 可视化结果
function visualize_results(nodes, displacements, strains, stresses)
    u = displacements(1:2:end);
    v = displacements(2:2:end);
    
    % 绘制位移
    figure;
    quiver(nodes(:, 1), nodes(:, 2), u, v);
    title('位移分布图');
    xlabel('x');
    ylabel('y');
    axis equal;
    
    % 绘制应变
    figure;
    contourf(nodes(:, 1), nodes(:, 2), reshape(strains(:, 1), size(nodes, 1), 1));
    colorbar;
    title('应变分布图');
    xlabel('x');
    ylabel('y');
    axis equal;
    
    % 绘制应力
    figure;
    contourf(nodes(:, 1), nodes(:, 2), reshape(stresses(:, 1), size(nodes, 1), 1));
    colorbar;
    title('应力分布图');
    xlabel('x');
    ylabel('y');
    axis equal;
end