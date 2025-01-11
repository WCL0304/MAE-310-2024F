% 组装单元刚度矩阵到全局刚度矩阵
function K = add_element_stiffness_to_global(K, nodes, k_elem)
    for i = 1:6
        for j = 1:6
            if k_elem(i, j) ~= 0
                K(nodes((i-1)*2 + 1), nodes((j-1)*2 + 1)) = K(nodes((i-1)*2 + 1), nodes((j-1)*2 + 1)) + k_elem(i, j);
                K(nodes((i-1)*2 + 2), nodes((j-1)*2 + 2)) = K(nodes((i-1)*2 + 2), nodes((j-1)*2 + 2)) + k_elem(i, j);
            end
        end
    end
end