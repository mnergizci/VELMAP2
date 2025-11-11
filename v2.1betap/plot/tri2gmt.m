function export_tri2gmt(matfile, outfile)
    % Export a triangular mesh to a GMT-readable multisegment file.
    % matfile: input .mat file containing the mesh (with a 'trim' struct)
    % outfile: output GMT multisegment file (e.g. 'triangles.xy')
    % Fengnian Chang

    % Load data
    data = load(matfile);
    if isfield(data, 'trim')
        tri = data.trim.tri;   % triangle connectivity (N x 3)
        x   = data.trim.x;     % x-coordinates of nodes
        y   = data.trim.y;     % y-coordinates of nodes
    else
        error('The input file does not contain the ''trim'' structure.');
    end

    % Write file
    fid = fopen(outfile, 'w');
    for i = 1:size(tri, 1)
        fprintf(fid, ">\n"); % each triangle as one segment
        for j = 1:3
            fprintf(fid, "%f %f\n", x(tri(i, j)), y(tri(i, j)));
        end
        % close the triangle
        fprintf(fid, "%f %f\n", x(tri(i, 1)), y(tri(i, 1)));
    end
    fclose(fid);

    fprintf('Exported %d triangles to %s\n', size(tri,1), outfile);
end
