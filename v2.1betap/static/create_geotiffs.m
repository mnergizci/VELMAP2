function create_geotiffs(insarfit, insarpar, sboifit, sboipar, outdirfile)
% Export GeoTIFFs from existing insarfit/sboifit preserving their (downsampled) size
% by building a new referencing object over the same geographic extent as vstd.


%% InSAR
insar_out_dir = fullfile(outdirfile, 'geotiffs');
if ~exist(insar_out_dir, 'dir'); mkdir(insar_out_dir); end

if ~isempty(insarfit)
    for i = 1:numel(insarfit)
        frame_id = get_frame_id(insarpar.dir{i});
        vstd_path = fullfile(insarpar.dir{i}, [frame_id, '.vstd_scaled.geo.tif']);
        if ~isfile(vstd_path)
            warning('Missing vstd: %s — skip InSAR %d', vstd_path, i); continue;
        end

        [vstd, Rvstd] = geotiffread(vstd_path);
        fit = insarfit(i);

        % Downsampled products
        grid    = -1 * (fit.ratemap + fit.resmap);
        model   = -1 *  fit.ratemap;
        resmap  = -1 *  fit.resmap;
        orbmap  = getfield_safe(fit,'orbmap',[]);
        atmmap  = getfield_safe(fit,'atmmap',[]);

        % Build new ref spanning same extent as vstd, but with size(grid)
        Rout = georefcells(Rvstd.LatitudeLimits, Rvstd.LongitudeLimits, size(grid), ...
                           'ColumnsStartFrom', Rvstd.ColumnsStartFrom, ...
                           'RowsStartFrom',    Rvstd.RowsStartFrom);

        % Try to copy CRS tags from vstd
        info = try_geotiffinfo(vstd_path);

        write_gt(insar_out_dir, [frame_id, '.vel_filt.mskd.eurasia.geo.tif'], grid,  Rout, info);
        % write_gt(insar_out_dir, [frame_id, '.vstd_scaled.geo.tif'],              vstd,  Rvstd, info);
        % if ~isempty(orbmap), write_gt(insar_out_dir, [frame_id, '_orbmap.tif'],  orbmap, Rout, info); end
        % if ~isempty(atmmap), write_gt(insar_out_dir, [frame_id, '_atmmap.tif'],  atmmap, Rout, info); end
        % write_gt(insar_out_dir, [frame_id, '_resmap.tif'],   resmap, Rout, info);
        write_gt(insar_out_dir, [frame_id, '_modelmap.tif'], model,  Rout, info);
    end
end

%% SBOI
sboi_out_dir = fullfile(outdirfile, 'geotiffs-sboi');
if ~exist(sboi_out_dir, 'dir'); mkdir(sboi_out_dir); end

if ~isempty(sboifit)
    for i = 1:numel(sboifit)
        frame_id = get_frame_id(sboipar.dir{i});
        vstd_path = fullfile(sboipar.dir{i}, [frame_id, '.vstd_scaled.geo.tif']);
        if ~isfile(vstd_path)
            warning('Missing vstd: %s — skip SBOI %d', vstd_path, i); continue;
        end

        [vstd, Rvstd] = geotiffread(vstd_path);
        fit = sboifit(i);

        grid    = -1 * (fit.ratemap + fit.resmap);
        model   = -1 *  fit.ratemap;
        resmap  = -1 *  fit.resmap;
        orbmap  = getfield_safe(fit,'orbmap',[]);
        atmmap  = getfield_safe(fit,'atmmap',[]);

        Rout = georefcells(Rvstd.LatitudeLimits, Rvstd.LongitudeLimits, size(grid), ...
                           'ColumnsStartFrom', Rvstd.ColumnsStartFrom, ...
                           'RowsStartFrom',    Rvstd.RowsStartFrom);

        info = try_geotiffinfo(vstd_path);

        write_gt(sboi_out_dir, [frame_id, '.vel_filt.mskd.eurasia.geo.tif'], grid,  Rout, info);
        % write_gt(sboi_out_dir, [frame_id, '.vstd_scaled.geo.tif'],              vstd,  Rvstd, info);
        % if ~isempty(orbmap), write_gt(sboi_out_dir, [frame_id, '_orbmap.tif'],  orbmap, Rout, info); end
        % if ~isempty(atmmap), write_gt(sboi_out_dir, [frame_id, '_atmmap.tif'],  atmmap, Rout, info); end
        % write_gt(sboi_out_dir, [frame_id, '_resmap.tif'],   resmap, Rout, info);
        write_gt(sboi_out_dir, [frame_id, '_modelmap.tif'], model,  Rout, info);
    end
end

%% ---------- helpers ----------
function v = getfield_safe(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function frame_id = get_frame_id(dirpath)
    parts = strsplit(dirpath, filesep);
    frame_id = parts{end-1};
end

function info = try_geotiffinfo(path)
    try
        info = geotiffinfo(path);
    catch
        info = [];
    end
end

function write_gt(folder, fname, A, R, info)
    if ~exist(folder,'dir'), mkdir(folder); end
    outpath = fullfile(folder, fname);
    if ~isempty(info) && isfield(info,'GeoTIFFTags') && isfield(info.GeoTIFFTags,'GeoKeyDirectoryTag')
        geotiffwrite(outpath, A, R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
    else
        geotiffwrite(outpath, A, R);
    end
end

end
