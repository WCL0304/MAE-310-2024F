% 计算应变
function strains = calculate_strains(displacements, nodes, elements)
    nElements = size(elements, 1);
    strains = zeros(nElements, 3);
    for i = 1:nElements
        element = elements(i, :);
        disp_elem = displacements(2 * (element - 1) + [1, 2], :);
        [B, detJ] = strain_displacement_matrix(element, nodes, [0, 0]);
        strains(i, :) = B * disp_elem / detJ;
    end
end