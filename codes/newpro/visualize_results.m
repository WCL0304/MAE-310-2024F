function visualize_results(nodes, displacements, strains, stresses)
    % 绘制位移
    figure;
    quiver(nodes(:, 1), nodes(:, 2), displacements(:, 1), displacements(:, 2));
    title('位移分布图');
    xlabel('x');
    ylabel('y');
    axis equal;

    % 绘制应变
    figure;
    contourf(nodes(:, 1), nodes(:, 2), strains);
    colorbar;
    title('应变分布图');
    xlabel('x');
    ylabel('y');
    axis equal;

    % 绘制应力
    figure;
    contourf(nodes(:, 1), nodes(:, 2), stresses);
    colorbar;
    title('应力分布图');
    xlabel('x');
    ylabel('y');
    axis equal;
end