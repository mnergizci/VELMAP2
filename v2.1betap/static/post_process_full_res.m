%=================================================================
% post_process_full_res.m
% 
% Muhammet Nergizci @ Leeds, 30/05/2025
%=================================================================

function post_process_full_res(outdirfile)
tic
%% Load inputs

% Construct full path to the file
fullpath = fullfile(outdirfile, 'precrash.mat');

% Check and load
if exist(fullpath, 'file') == 2
    vars = load(fullpath);
    % Unpack variables
    trim     = vars.trim;
    fitmodel = vars.fitmodel;
    gps      = vars.gps;
    invenu   = vars.invenu;
    outdir   = vars.outdir;
    insarpar = vars.insarpar;
    sboipar  = vars.sboipar;
else
    error('Cannot find file: %s\nCheck if the path is correct or if the file exists.', fullpath);
end

%% 1x1 looking forward calculation for InSAR
if insarpar.ninsarfile > 0 && insarpar.activate > 0
    insar_1lk_file = fullfile(outdirfile, 'insar_1lk.mat');
    insarfit1lk_file = fullfile(outdirfile, 'insarfit1lk.mat');

    if exist(insar_1lk_file, 'file') ~= 2
        insar = loadlics(insarpar, 1); % Full-res load
        save(insar_1lk_file, 'insar', '-v7.3');
    else
        tmp = load(insar_1lk_file);
        insar = tmp.insar;
    end

    if exist(insarfit1lk_file, 'file') ~= 2
        insarfit1lk = insarfwd(insar, trim, fitmodel, invenu, outdir, gps);
        save(insarfit1lk_file, 'insarfit1lk', '-v7.3');
    else
        fprintf('Skipping InSAR forward modeling, file exists: %s\n', insarfit1lk_file);
        tmp = load(insarfit1lk_file);
        insarfit1lk = tmp.insarfit1lk;
    end
end

%% 1x1 looking forward calculation for SBOI
if sboipar.nsboifile > 0 && sboipar.activate > 0
    sboi_1lk_file = fullfile(outdirfile, 'sboi_1lk.mat');
    sboifit1lk_file = fullfile(outdirfile, 'sboifit1lk.mat');

    if exist(sboi_1lk_file, 'file') ~= 2
        sboi = loadlics_sboi(sboipar, 1); % Full-res load
        save(sboi_1lk_file, 'sboi', '-v7.3');
    else
        tmp = load(sboi_1lk_file);
        sboi = tmp.sboi;
    end

    if exist(sboifit1lk_file, 'file') ~= 2
        sboifit1lk = insarfwd(sboi, trim, fitmodel, invenu, outdir, gps, 1);
        save(sboifit1lk_file, 'sboifit1lk', '-v7.3');
    else
        fprintf('Skipping SBOI forward modeling, file exists: %s\n', sboifit1lk_file);
        tmp = load(sboifit1lk_file);
        sboifit1lk = tmp.sboifit1lk;
    end
end




%% Export full-resolution GeoTIFFs for InSAR
insar_out_dir = fullfile(outdirfile, 'geotiffs-1.1-m');
if ~exist(insar_out_dir, 'dir'); mkdir(insar_out_dir); end

if exist('insarfit1lk', 'var')
    for i = 1:length(insarfit1lk)
        % Extract frame name from insarpar.dir
        parts = strsplit(insarpar.dir{i}, filesep);
        frame_id = parts{end-1};

        % Path to vstd
%         vstd_path = fullfile(insarpar.dir{i}, 'vstd.geo.tif');
        vstd_path = fullfile(insarpar.dir{i}, [frame_id, '.vstd_scaled.geo.tif']);
        if ~isfile(vstd_path)
            warning('Missing vstd file: %s — skipping InSAR frame %d', vstd_path, i);
            continue;
        end

        [vstd, R] = geotiffread(vstd_path);
        fit = insarfit1lk(i);

        % Write outputs using frame_id in filenames
        write_map(insar_out_dir, [frame_id, '.vel_filt.mskd.eurasia.geo.tif'], -1 * (fit.ratemap + fit.resmap), R);
        % write_map(insar_out_dir, [frame_id, '.vstd_scaled.geo.tif'], vstd, R);
        % if isfield(fit, 'orbmap') && ~isempty(fit.orbmap)
        %     write_map(insar_out_dir, [frame_id, '_orbmap.tif'], fit.orbmap, R);
        % end
        % if isfield(fit, 'atmmap') && ~isempty(fit.atmmap)
        %     write_map(insar_out_dir, [frame_id, '_atmmap.tif'], fit.atmmap, R);
        % end
        % write_map(insar_out_dir, [frame_id, '_resmap.tif'], -1 * fit.resmap, R);
        write_map(insar_out_dir, [frame_id, '_modelmap.tif'], -1 * fit.ratemap, R);
    end
end

%% Export full-resolution GeoTIFFs for SBOI
sboi_out_dir = fullfile(outdirfile, 'geotiffs-1.1-m-sboi');
if ~exist(sboi_out_dir, 'dir'); mkdir(sboi_out_dir); end

if exist('sboifit1lk', 'var')
    for i = 1:length(sboifit1lk)
        parts = strsplit(sboipar.dir{i}, filesep);
        frame_id = parts{end-1};

        vstd_path = fullfile(sboipar.dir{i}, [frame_id, '.vstd_scaled.geo.tif']);
        if ~isfile(vstd_path)
            warning('Missing vstd file: %s — skipping SBOI frame %d', vstd_path, i); continue;
        end

        [vstd, R] = geotiffread(vstd_path);
        fit = sboifit1lk(i);

        grid = -1 * (fit.ratemap + fit.resmap);
        [ngrid_rows, ngrid_cols] = size(grid);
        [vstd_rows, vstd_cols] = size(vstd);

        if ~isequal([ngrid_rows, ngrid_cols], R.RasterSize)
            warning('Adjusting frame %s: size(grid) = [%d %d], R.RasterSize = [%d %d]', ...
                frame_id, ngrid_rows, ngrid_cols, R.RasterSize(1), R.RasterSize(2));

            % Use min size across all arrays to avoid index errors
            nrows = min([ngrid_rows, R.RasterSize(1), vstd_rows]);
            ncols = min([ngrid_cols, R.RasterSize(2), vstd_cols]);

            % Resize all arrays accordingly
            grid = grid(1:nrows, 1:ncols);
            vstd = vstd(1:nrows, 1:ncols);
            R.RasterSize = [nrows, ncols];

            if isfield(fit, 'ratemap')
                fit.ratemap = fit.ratemap(1:nrows, 1:ncols);
            end
            if isfield(fit, 'resmap')
                fit.resmap = fit.resmap(1:nrows, 1:ncols);
            end
            if isfield(fit, 'orbmap') && ~isempty(fit.orbmap)
                fit.orbmap = fit.orbmap(1:nrows, 1:ncols);
            end
            if isfield(fit, 'atmmap') && ~isempty(fit.atmmap)
                fit.atmmap = fit.atmmap(1:nrows, 1:ncols);
            end
        end

        write_map(sboi_out_dir, [frame_id, '.vel_filt.mskd.eurasia.geo.tif'], grid, R);
        % write_map(sboi_out_dir, [frame_id, '.vstd_scaled.geo.tif'], vstd, R);
        % 
        % if isfield(fit, 'orbmap') && ~isempty(fit.orbmap)
        %     write_map(sboi_out_dir, [frame_id, '_orbmap.tif'], fit.orbmap, R);
        % end
        % if isfield(fit, 'atmmap') && ~isempty(fit.atmmap)
        %     write_map(sboi_out_dir, [frame_id, '_atmmap.tif'], fit.atmmap, R);
        % end
        % 
        % write_map(sboi_out_dir, [frame_id, '_resmap.tif'], -1 * fit.resmap, R);
        write_map(sboi_out_dir, [frame_id, '_modelmap.tif'], -1 * fit.ratemap, R);
    end
end

function write_map(folder, filename, grid, R)
    outpath = fullfile(folder, filename);
    geotiffwrite(outpath, grid, R);
end

toc
end

