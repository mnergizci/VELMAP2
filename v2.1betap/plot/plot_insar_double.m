%=================================================================
% plot_insar.m
% 
% Tiledlayout requires 2019b or later.
%
% Modified from cahb_insarfit_milan.m.
% Andrew Watson @ leeds, 16/07/2021
% Qi Ou @ leeds, 20/06/2022
% Muhammet Nergizi @ leeds 26/04/2025
%=================================================================

function plot_insar_double(insarfit,gps,outdir, pngfile)

if nargin < 4 || isempty(pngfile)
    pngfile = 'model_los.png';
end

% parameters
clim = [-20 20];
places = [];%{'Iran Islamic Republic of','Iraq','Afghanistan','Turkey'}
bordersfile = []; %'borderdata.mat';

%% load inputs

% plotting
S2 = load('vik.mat');     % S2 is just an example name
vik = S2.vik;             % now assign the colormap
S3= load('borderdata.mat');   %
borders = S3.places;

%% format inputs
% ---- pre-alloc & fill
n = numel(insarfit);
x = cell(1,n); y = cell(1,n);
datamap = cell(1,n); refmap = cell(1,n);
resimap = cell(1,n); ratemap = cell(1,n);
orbmap  = cell(1,n); atmmap  = cell(1,n);
insarboundary = cell(1,n); passdir = cell(1,n);



for ii = 1:n
    % coords
    x{ii} = insarfit(ii).ifghdr.xfirst + (0:insarfit(ii).ifghdr.width-1).*insarfit(ii).ifghdr.xstep;
    y{ii} = insarfit(ii).ifghdr.yfirst + (0:insarfit(ii).ifghdr.length-1).*insarfit(ii).ifghdr.ystep;

    % frame boundary (on valid pixels only)
    [xx, yy] = meshgrid(x{ii}, y{ii});
    xx = xx(:); yy = yy(:);
    mask = ~isnan(insarfit(ii).ratemap(:));
    xx = xx(mask); yy = yy(mask);
    if numel(xx) >= 3
        boundaryind = boundary(xx(:), yy(:));
        insarboundary{ii} = [xx(boundaryind) yy(boundaryind)];
    else
        insarboundary{ii} = [NaN NaN];
    end

    passdir{ii} = insarfit(ii).ifghdr.passdir;

    % maps
    datamap{ii} = insarfit(ii).stackmap + insarfit(ii).resmap; % observed LOS
    refmap{ii}  = insarfit(ii).ratemap + insarfit(ii).resmap;  % referenced LOS
    resimap{ii} = insarfit(ii).resmap;
    stackmap{ii}= insarfit(ii).stackmap;
    ratemap{ii} = insarfit(ii).ratemap;                         % model LOS
    orbmap{ii}  = insarfit(ii).orbmap;                          % model ramps
    atmmap{ii}  = insarfit(ii).stackmap - insarfit(ii).orbmap - insarfit(ii).ratemap; % model atm
end

% ----------------------------------------------------------------
% Asc/Desc indices
% ----------------------------------------------------------------
asc_ind  = strcmp(passdir, 'A');
desc_ind = strcmp(passdir, 'D');

% ----------------------------------------------------------------
% Color limits (robust percentiles per pass)
% ----------------------------------------------------------------
[clim_asc, clim_desc] = get_clims(refmap, asc_ind, desc_ind);

% % If original default output name, force symmetric ±20 mm/yr like before
% if strcmp(pngfile, "model_los.png")
%     clim_asc  = [-20 20];
%     clim_desc = [-20 20];
% end

%% colorbar limits
ratemap_asc = [];
ratemap_desc = [];
for ii = 1:length(ratemap)
    if asc_ind(ii)
        ratemap_asc = [ratemap_asc; ratemap{ii}(:)];
    elseif desc_ind(ii)
        ratemap_desc = [ratemap_desc; ratemap{ii}(:)];
    end
end

ratemap_asc = ratemap_asc(~isnan(ratemap_asc));
ratemap_desc = ratemap_desc(~isnan(ratemap_desc));

cmin_asc = round(prctile(ratemap_asc,1));
cmax_asc = round(prctile(ratemap_asc,99));
cmin_desc = round(prctile(ratemap_desc,1));
cmax_desc = round(prctile(ratemap_desc,99));

clim_asc = [cmin_asc cmax_asc];
clim_desc = [cmin_desc cmax_desc];
% 
% if pngfile == "model_los.png"
%     clim_asc = [-100 100];
%     clim_desc = [-100 100];
% end 

% ----------------------------------------------------------------
% Axes limits
% ----------------------------------------------------------------
lonlim = [min(cell2mat(x)), max(cell2mat(x))];
latlim = [min(cell2mat(y)), max(cell2mat(y))];

% ----------------------------------------------------------------
% Make two figures: Ascending and Descending
% ----------------------------------------------------------------
[~, base, ~] = fileparts(pngfile);
png_asc = fullfile(outdir, [base '_asc.png']);
png_dsc = fullfile(outdir, [base '_dsc.png']);

make_panel_figure(x, y, refmap, datamap, ratemap, orbmap, resimap, atmmap, ...
    insarboundary, asc_ind, clim_asc, borders, places, vik, lonlim, latlim, ...
    ['Ascending — ' outdir], png_asc);

make_panel_figure(x, y, refmap, datamap, ratemap, orbmap, resimap, atmmap, ...
    insarboundary, desc_ind, clim_desc, borders, places, vik, lonlim, latlim, ...
    ['Descending — ' outdir], png_dsc);

end % ====== main function ======

%=================================================================
% Helpers
%=================================================================

function [clim_asc, clim_desc] = get_clims(refmap, asc_ind, desc_ind)
    refmap_asc  = [];
    refmap_desc = [];
    for ii = 1:length(refmap)
        if asc_ind(ii)
            refmap_asc  = [refmap_asc;  refmap{ii}(:)];
        elseif desc_ind(ii)
            refmap_desc = [refmap_desc; refmap{ii}(:)];
        end
    end
    refmap_asc  = refmap_asc(~isnan(refmap_asc));
    refmap_desc = refmap_desc(~isnan(refmap_desc));

    if isempty(refmap_asc),  refmap_asc = 0;  end
    if isempty(refmap_desc), refmap_desc = 0; end

    cmin_asc  = round(prctile(refmap_asc,  1));
    cmax_asc  = round(prctile(refmap_asc, 99));
    cmin_desc = round(prctile(refmap_desc, 1));
    cmax_desc = round(prctile(refmap_desc,99));

    clim_asc  = [cmin_asc  cmax_asc];
    clim_desc = [cmin_desc cmax_desc];
end

function make_panel_figure(x, y, refmap, datamap, ratemap, orbmap, resimap, atmmap, ...
    insarboundary, pass_idx, clim_pass, borders, places, cmap, lonlim, latlim, ...
    title_pass, outfile)

    figure('Units','Normalized','OuterPosition',[0.2,0.2,0.6,0.55]);
    tiledlayout(2,3,'TileSpacing','compact');

    % 1) Referenced LOS
    nexttile; hold on
    plot_data(x(pass_idx), y(pass_idx), refmap(pass_idx), insarboundary(pass_idx), ...
        borders, places, lonlim, latlim, clim_pass);
    title('Referenced LOS');

    % 2) Observed LOS (stack + res)
    nexttile; hold on
    plot_data(x(pass_idx), y(pass_idx), datamap(pass_idx), insarboundary(pass_idx), ...
        borders, places, lonlim, latlim, clim_pass);
    title('Observed LOS');

    % 3) Model LOS (rate map)
    nexttile; hold on
    plot_data(x(pass_idx), y(pass_idx), ratemap(pass_idx), insarboundary(pass_idx), ...
        borders, places, lonlim, latlim, clim_pass);
    title('Model LOS');

    % 4) Model ramps (orbital)
    nexttile; hold on
    plot_data(x(pass_idx), y(pass_idx), orbmap(pass_idx), insarboundary(pass_idx), ...
        borders, places, lonlim, latlim, clim_pass);
    title('Model Ramps');

    % 5) Residual
    nexttile; hold on
    plot_data(x(pass_idx), y(pass_idx), resimap(pass_idx), insarboundary(pass_idx), ...
        borders, places, lonlim, latlim, clim_pass);
    title('Residual');

    % 6) Modeled atmosphere
    nexttile; hold on
    plot_data(x(pass_idx), y(pass_idx), atmmap(pass_idx), insarboundary(pass_idx), ...
        borders, places, lonlim, latlim, clim_pass);
    title('Modeled Atmosphere');

    colormap(cmap)
    sgtitle(strrep(title_pass, '_', ' '));

    % One colorbar overall (attach to last axes)
    cb = colorbar; cb.Label.String = 'mm/yr';

    saveas(gcf, outfile);
end

function plot_data(xc, yc, data, insarboundary, borders, places, lonlim, latlim, clim)
    for jj = 1:length(xc)
        imagesc(xc{jj}, yc{jj}, data{jj}, 'AlphaData', ~isnan(data{jj}));
    end

    % Optional borders
    if ~isempty(borders) && ~isempty(places)
        for ii = 1:length(places)
            b_ind = find(strcmp(borders.places, places(ii)));
            if ~isempty(b_ind)
                plot(borders.lon{b_ind}, borders.lat{b_ind}, 'r');
            end
        end
    end

    set(gca,'YDir','normal'); % ensure north-up
    caxis(clim); xlim(lonlim); ylim(latlim);
    axis tight; axis equal;
end
