function trim = mesh4velmap(outlinefile, dx, dy, outfile)
% mesh4velmap  Create a smoothed triangular mesh inside a polygon
% 
%   trim = mesh4velmap(outlinefile, dx, dy)
%   trim = mesh4velmap(outlinefile, dx, dy, outfile)
%
% Inputs
%   outlinefile : text file with polygon vertices [x y] (e.g. 'mesh_poly.gmt')
%   dx, dy      : grid spacing in x and y (e.g. 0.2, 0.2)
%   outfile     : optional .mat file to save mesh (default: 'mesh_<dx>.mat')
%
% Output
%   trim : struct with fields
%          trim.x   - node x coordinates
%          trim.y   - node y coordinates
%          trim.tri - triangle connectivity (Nx3)
%
% Jin Fang, COMET-University of Leeds, 24/11/2025

    %------------------------------
    % defaults and setup
    %------------------------------
    if nargin < 3
        error('Usage: mesh4velmap(outlinefile, dx, dy, [outfile])');
    end

    if nargin < 4 || isempty(outfile)
        % default name based on dx
        outfile = sprintf('mesh_%g.mat', dx);
    end

    meshoutline = 1;   % always use outline
    alpha_val   = 1;   % alpha shape parameter

    % add paths if needed (adjust to your directory structure)
    addpath ./mesh2d
    addpath ./tools

    %------------------------------
    % read outline polygon
    %------------------------------
    outline = load(outlinefile);   % expects [x y] columns

    xmin = min(outline(:,1)) - 0.5;
    xmax = max(outline(:,1)) + 0.5;
    ymin = min(outline(:,2)) - 0.5;
    ymax = max(outline(:,2)) + 0.5;

    %------------------------------
    % generate regular grid of nodes
    %------------------------------
    [x, y] = meshgrid(xmin:dx:xmax, ymin:dy:ymax);
    M_nodes = [x(:) y(:)];

    M_nodes_X = M_nodes(:,1);
    M_nodes_Y = M_nodes(:,2);

    %------------------------------
    % keep only nodes inside polygon
    %------------------------------
    if meshoutline == 1
        [in_poly, ~] = inpolygon(M_nodes_X, M_nodes_Y, outline(:,1), outline(:,2));
        M_nodes_X(~in_poly) = [];
        M_nodes_Y(~in_poly) = [];
        M_nodes = [M_nodes_X M_nodes_Y];
    end

    %------------------------------
    % alpha-shape triangulation + smoothing
    %------------------------------
    shp = alphaShape(M_nodes(:,1), M_nodes(:,2), alpha_val);
    tri = alphaTriangulation(shp);

    [edge] = tricon2(tri);
    ebnd = edge(:,4) < +1;        % boundary edges
    conn = edge(ebnd, 1:2);

    [M_nodes, etri, trim0_tri, tnum] = smooth2(M_nodes, conn, tri); %#ok<ASGLU>

    %------------------------------
    % build trim structure
    %------------------------------
    trim0.x   = M_nodes(:,1);
    trim0.y   = M_nodes(:,2);
    trim0.tri = trim0_tri;

    % preview
    figure;
    patch('faces', trim0.tri(:,1:3), 'vertices', M_nodes, ...
          'facecolor', 'w', ...
          'edgecolor', [.6 .6 .6]);
    axis equal
    title('mesh4velmap mesh')

    saveas(gcf, 'MESH.png');

    trim = trim0;
    save(outfile, 'trim');
end
