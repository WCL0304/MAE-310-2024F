% 计算应力
function stresses = calculate_stresses(strains, E, nu, isPlaneStress)
    D = plane_stress_strain_matrix(E, nu, isPlaneStress);
    stresses = D * strains';
end