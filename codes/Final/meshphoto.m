% % 加载网格数据
% try
%     load('mesh.m'); % 确保文件在当前目录下
% catch ME
%     warning('数据文件 "mesh.m" 未找到，请确保文件存在于当前目录。');
%     return;
% end

% 提取节点和单元信息
nodeCoordinates = msh.POS(:, 1:2); % 节点坐标 (x, y)
triangleElements = msh.TRIANGLES(:, 1:3); % 三角形单元 (节点索引)

% 绘制网格
figure;
trimesh(triangleElements, nodeCoordinates(:, 1), nodeCoordinates(:, 2), 'EdgeColor', 'k');
axis equal;
title('四分之一带孔平板网格');
xlabel('X');
ylabel('Y');