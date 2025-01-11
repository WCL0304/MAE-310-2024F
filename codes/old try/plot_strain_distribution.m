% 绘制应变分布到网格
function plot_strain_distribution(IEN, x_coor, y_coor, strain_x)
    figure;
    patch('Faces', IEN, 'Vertices', [x_coor, y_coor], 'FaceVertexCData', strain_x, ...
        'FaceColor', 'flat', 'EdgeColor', 'k');
    colorbar;
    title('四分之一带孔平板的应变分布 (εx)');
    xlabel('X 坐标');
    ylabel('Y 坐标');
    axis equal;
end