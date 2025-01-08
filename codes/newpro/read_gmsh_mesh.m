function [nodes, elements, boundaryEdges] = read_gmsh_mesh(filename)
    % 读取Gmsh网格数据
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件：%s', filename);
    end

    % 读取节点数据
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line) && startsWith(line, '$Nodes')
            nNodes = fscanf(fid, '%d', 1);
            disp(['读取的节点数量: ', num2str(nNodes)]); % 调试信息
            if isnan(nNodes) || isinf(nNodes)
                error('节点数量读取错误，值为 NaN 或 Inf');
            end
            nodes = zeros(nNodes, 2);
            for i = 1:nNodes
                data = fscanf(fid, '%d %f %f %f', [4, 1]);
                nodes(i, :) = data(2:3);
            end
            break;
        end
    end

    % 读取元素数据
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line) && startsWith(line, '$Elements')
            nElements = fscanf(fid, '%d', 1);
            disp(['读取的元素数量: ', num2str(nElements)]); % 调试信息
            elements = [];
            boundaryEdges = [];
            for i = 1:nElements
                data = fscanf(fid, '%d %d %d %d %d %d', [6, 1]); % 读取6个整数
                if data(2) == 2 % 三角形单元
                    elements = [elements; data(4:6)];
                elseif data(2) == 15 % 边界线
                    if length(data) < 6
                        error('边界元素数据格式错误，数据长度不足');
                    end
                    boundaryEdges = [boundaryEdges; data(4:5)];
                end
            end
            break;
        end
    end

    % 检查 boundaryEdges 的格式
    if size(boundaryEdges, 2) ~= 2
        error('边界元素数组格式错误，每行应包含两个节点索引');
    end

    fclose(fid);
end