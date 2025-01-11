% 计算单元刚度矩阵
function K_elem = calculate_element_stiffness_matrix(element, nodes, E, nu, isPlaneStress)
    [B, detJ] = strain_displacement_matrix(element, nodes, [0, 0]);
    D = plane_stress_strain_matrix(E, nu, isPlaneStress);
    K_elem = detJ * B' * D * B;
end