% 创建 LM 数组 (局部到全局自由度映射)
function LM = create_LM(IEN, ID)
    n_elements = size(IEN, 1);
    LM = zeros(n_elements, 6); % 3 节点 * 2 DOFs/节点
    for ee = 1:n_elements
        for aa = 1:3
            node = IEN(ee, aa);
            LM(ee, (aa-1)*2 + 1) = ID(node, 1);
            LM(ee, (aa-1)*2 + 2) = ID(node, 2);
        end
    end
end