% 计算应变-位移矩阵和雅各比矩阵的行列式
function [B, detJ] = strain_displacement_matrix(element, nodes, gauss_point)
    % 计算应变-位移矩阵和雅各比矩阵的行列式
    % 假设使用线性三角形单元
    p1 = nodes(element(1), :);
    p2 = nodes(element(2), :);
    p3 = nodes(element(3), :);
    
    % 雅各比矩阵
    J = [p2(1) - p1(1), p3(1) - p1(1);
         p2(2) - p1(2), p3(2) - p1(2)];
    detJ = det(J);
    
    % 应变-位移矩阵
    B = [p2(2) - p3(2), 0, p3(2) - p1(2), 0, p1(2) - p2(2), 0;
         0, p3(1) - p2(1), 0, p1(1) - p3(1), 0, p2(1) - p1(1);
         p2(2) - p3(2), p3(1) - p2(1), p3(2) - p1(2), p1(1) - p3(1), p1(2) - p2(2), p2(1) - p1(1)] / (2 * detJ);
end