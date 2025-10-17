function exportToVTK(vertices, faces, filename, varargin)
% exportToVTK Exports vertices and optionally faces to a VTK file in legacy format.
%
%   exportToVTK(vertices, faces, filename)
%   exportToVTK(vertices, faces, filename, 'PointsOnly', true)
%
%   vertices: an N x 3 matrix of vertex coordinates.
%   faces: an M x K matrix of face indices (typically 3 for triangles).
%          Can be empty if not applicable.
%   filename: the output file name (e.g., 'mesh.vtk').
%
%   Optional parameter:
%       'PointsOnly' (default false): if true, only the points are saved.
%           The points are written as individual vertex cells.
%
% Example:
%   % Export mesh (points and faces):
%   vertices = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
%   faces = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
%   exportToVTK(vertices, faces, 'mesh.vtk');
%
%   % Export only points:
%   exportToVTK(vertices, faces, 'pointsOnly.vtk', 'PointsOnly', true);

    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'PointsOnly', false, @islogical);
    parse(p, varargin{:});
    pointsOnly = p.Results.PointsOnly;
    
    fid = fopen(filename, 'w');
    if fid < 0
        error('Could not open file %s for writing.', filename);
    end

    % Write VTK header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'VTK output from MATLAB\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET POLYDATA\n');

    % Write the points section
    numPoints = size(vertices, 1);
    fprintf(fid, 'POINTS %d float\n', numPoints);
    for i = 1:numPoints
        fprintf(fid, '%f %f %f\n', vertices(i,1), vertices(i,2), vertices(i,3));
    end
    
    if pointsOnly
        % Write the vertices as individual vertex cells
        % Each vertex cell contains 2 numbers: the count (1) and the vertex index.
        totalVertexData = numPoints * 2;
        fprintf(fid, 'VERTICES %d %d\n', numPoints, totalVertexData);
        for i = 1:numPoints
            fprintf(fid, '1 %d\n', i - 1);  % Convert MATLAB indices (1-indexed) to VTK (0-indexed)
        end
    else
        % If faces are provided, write them as polygons
        if ~isempty(faces)
            numFaces = size(faces, 1);
            verticesPerFace = size(faces, 2);
            totalFaceData = numFaces * (verticesPerFace + 1);
            fprintf(fid, 'POLYGONS %d %d\n', numFaces, totalFaceData);
            for i = 1:numFaces
                fprintf(fid, '%d ', verticesPerFace);
                fprintf(fid, '%d ', faces(i,:) - 1);  % Convert to 0-indexed
                fprintf(fid, '\n');
            end
        end
    end

    fclose(fid);
    fprintf('VTK file saved as %s\n', filename);
end
