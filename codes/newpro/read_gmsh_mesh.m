function [nodes, elements, boundaryEdges] = read_gmsh_mesh(filename)
    % 读取 Gmsh 网格数据
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件: %s', filename);
    end

    % 读取节点数据
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line) && startsWith(line, '$Nodes')
            nNodes = fscanf(fid, '%d', 1);
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
            elements = [];
            boundaryEdges = [];
            for i = 1:nElements
                data = fscanf(fid, '%d %d %d %d %d %d %d', [7, 1]);
                if data(2) == 2 % 三角形单元
                    elements = [elements; data(5:7)];
                elseif data(2) == 15 % 边界线
                    boundaryEdges = [boundaryEdges; data(5:6)];
                end
            end
            break;
        end
    end

    fclose(fid);
end