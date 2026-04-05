function [xyz, faces] = read_ply_xyz(filename)
% Reads a PLY file and returns Nx3 points (vertices) and optional faces
fid = fopen(filename, 'r');
if fid == -1, error('Cannot open file: %s', filename); end

% Parse header
nverts = 0; nfaces = 0; headerLines = 0;
while true
    tline = fgetl(fid);
    headerLines = headerLines + 1;
    if contains(tline, 'element vertex')
        nverts = sscanf(tline, 'element vertex %d');
    elseif contains(tline, 'element face')
        nfaces = sscanf(tline, 'element face %d');
    elseif strcmp(tline, 'end_header')
        break;
    end
end

% Read vertices
xyz = fscanf(fid, '%f %f %f', [3, nverts])';
faces = [];

% If there are faces, read them (optional)
if nfaces > 0
    fseek(fid, 0, 'bof'); % rewind
    % skip header + vertex lines
    for i = 1:headerLines+nverts
        fgetl(fid);
    end
    faces = zeros(nfaces, 3);
    for i = 1:nfaces
        data = fscanf(fid, '%d', 1+3); % first value is number of vertices for the face (should be 3)
        if isempty(data) || numel(data) < 4
            continue;
        end
        faces(i,:) = data(2:4)' + 1; % convert to 1-based indexing
    end
end
fclose(fid);
end