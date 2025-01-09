% 读取Gmsh网格数据
function [nodes, elements, boundaryEdges] = read_gmsh_mesh(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件：%s', filename);
    end
    
    % 读取节点数据
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line) && startsWith(line, '$Nodes')
            nNodes = fscanf(fid, '%d', 1);
            disp(['读取的节点数量:', num2str(nNodes)]);
            if isnan(nNodes) || isinf(nNodes)
                error('节点数量读取错误，值为NaN或Inf');
            end
            nodes = zeros(nNodes, 2);
            for i = 1:nNodes
                node_data = fscanf(fid, '%d %f %f', [3, 1]);
                nodes(i, :) = node_data(2:3);
            end
            line = fgetl(fid); % 用于忽略节点数据后的结束符
            break;
        end
    end
    
    % 读取元素数据
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line) && startsWith(line, '$Elements')
            nElements = fscanf(fid, '%d', 1);
            disp(['读取的元素数量:', num2str(nElements)]);
            elements = [];
            boundaryEdges = [];
            for i = 1:nElements
                element_data = fscanf(fid, '%d %d %d %d %d %d', [6, 1]);
                if element_data(2) == 2 % 三角形单元
                    elements = [elements; element_data(4:6)];
                elseif element_data(2) == 15 % 边界线
                    if length(element_data) < 6
                        error('边界元素数据格式错误，数据长度不足');
                    end
                    boundaryEdges = [boundaryEdges; element_data(4:5)];
                end
            end
            line = fgetl(fid); % 用于忽略元素数据后的结束符
            break;
        end
    end
    
    % 检查boundaryEdges的格式
    if size(boundaryEdges, 2) ~= 2
        error('边界元素数组格式错误，每行应包含两个节点索引');
    end
    
    fclose(fid);
end