% 计算平面应力/应变矩阵
function D = plane_stress_strain_matrix(E, nu, isPlaneStress)
    % 平面应力/应变矩阵
    if isPlaneStress
        D = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    else
        D = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
    end
end