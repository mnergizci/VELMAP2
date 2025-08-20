%=================================================================
% myplot_vel_strain.m
% Plot velocities and strains from velmap output.
% Tiledlayout requires 2019b or later.
%
% Modified from cahb_insarfit_gpsonly.m.
% Andrew Watson @ leeds, 12/07/2021
% Qi Ou @ leeds, 23/06/2023 
% Jin Fang @ Leeds, 23/11/2024
% Muhammet Nergizci @ leeds, 26/04/2025
%=================================================================

function myplot_vel_strain(outdir, gps, outdir2)

if nargin < 2 || isempty(gps)
    if exist(fullfile(outdir, 'gps.mat'), 'file') == 2
        S = load(fullfile(outdir, 'gps.mat'));
        gps = S.gps;  % load variable inside
    else
        warning('gps.mat not found. GPS data will not be plotted.');
        gps = [];
    end
end

if nargin < 3 || isempty(outdir2)
    plot_diff = false;
else
    plot_diff = true;
end

% main paths
velfile = strcat(outdir,'velfit.dat');
strainfile = strcat(outdir,'strain_savage_nring1.dat');
faultfile = 'faults/gem_active_faults.gmt';

%second paths
if plot_diff
    velfile2 = strcat(outdir2,'velfit.dat');
    strainfile2 = strcat(outdir2,'strain_savage_nring1.dat');
end

% parameters
meshspacing = 0.01;
places = {'China', 'Kazakhstan', 'Kyrgyzstan', 'Uzbekistan', 'Tajikistan'}; % to plot

%% load inputs

% results
% gps = []; %readmatrix(gpspar.filename{1});
vel = readmatrix(velfile);
strain = readmatrix(strainfile);

%read second inputs
if plot_diff
    vel2 = readmatrix(velfile2);
    strain2 = readmatrix(strainfile2);
end

% plotting
Sborder= load('borderdata.mat');
borders = Sborder.places;
Svik = load('vik.mat');
vik = Svik.vik;
Svik = load('vik.mat');
vik = Svik.vik;
Sacton = load('acton.mat');
acton = flipud(Sacton.acton);
Slajolla = load('lajolla.mat');
lajolla = flipud(Slajolla.lajolla);
Sbam = load('bam.mat');
bam = Sbam.bam;

% plotting coords
lon = vel(:,1); lat = vel(:,2);
%poly = readmatrix(intersection_txt);
% faults
try
    load faults.mat;
catch
    faults = read_faults_from_gmt(faultfile)
    save faults.mat faults;
end


%% format velocities
%first enu
[egrid,ngrid,ugrid] = format_vel(vel,meshspacing);
%std
[estdgrid,nstdgrid,usdtgrid] = format_std(vel,meshspacing);

%second enu
if plot_diff
    [egrid2,ngrid2,ugrid2] = format_vel(vel2,meshspacing);
end

% JF
lon_grid = linspace(min(lon), max(lon), size(egrid, 2)); 
lat_grid = linspace(min(lat), max(lat), size(egrid, 1));

%% format strain
%first strain
[strainmap,digrid,ssgrid,msgrid,i2grid] = format_strain(strain,meshspacing);
%second strain
if plot_diff
    [strainmap2,digrid2,ssgrid2,msgrid2,i2grid2] = format_strain(strain2,meshspacing);
end


%% plot vels

clim_straim = [0 1e-7];

if ~plot_diff
    %%%%
    % East
    data_flat_e = egrid(:);
    data_flat_e = data_flat_e(~isnan(data_flat_e));
    clim_e = round([prctile(data_flat_e, 1) prctile(data_flat_e, 99)]);
    
    % North
    data_flat_n = ngrid(:);
    data_flat_n = data_flat_n(~isnan(data_flat_n));
    clim_n = round([prctile(data_flat_n, 1) prctile(data_flat_n, 99)]);
    
    % Up
    data_flat_u = ugrid(:);
    data_flat_u = data_flat_u(~isnan(data_flat_u));
    clim_u = round([prctile(data_flat_u, 1) prctile(data_flat_u, 99)]);
    %%%%
    
    
    figure()
    tiledlayout(2,3,'TileSpacing','compact')
    axis equal;
    
    ax = nexttile; hold on
    % plt_data(lon,lat,egrid,gps,strainmap.lonlims,strainmap.latlims,[-15 15],borders,places,'East vel (mm/yr)', faults, poly)
    %JF
    plt_data(lon_grid,lat_grid,egrid,gps,strainmap.lonlims,strainmap.latlims,clim_e,borders,places,'East vel (mm/yr)', faults, 0)%, poly)
    colormap(ax, vik)
    
    ax = nexttile; hold on
    % plt_data(lon,lat,ngrid,gps,strainmap.lonlims,strainmap.latlims,[-30, 30],borders,places,'North vel (mm/yr)', faults)%, poly)
    %JF
    plt_data(lon_grid,lat_grid,ngrid,gps,strainmap.lonlims,strainmap.latlims,clim_n,borders,places,'North vel (mm/yr)', faults, 0)%, poly)
    colormap(ax, vik)
    
    ax = nexttile; hold on
    % plt_data(lon,lat,ugrid,gps,strainmap.lonlims,strainmap.latlims,[-15, 15],borders,places,'Up vel (mm/yr)', faults)%, poly)
    %JF
    plt_data(lon_grid,lat_grid,ugrid,gps,strainmap.lonlims,strainmap.latlims,clim_u,borders,places,'Up vel (mm/yr)', faults, 1)%, poly)
    colormap(ax, vik)
    
    ax = nexttile; hold on
    % plt_data(lon,lat,msgrid,gps,strainmap.lonlims,strainmap.latlims,clim_straim,borders,places,...
    %     'Max shear (/yr)', faults)%, poly)
    %JF
    plt_data(lon_grid,lat_grid,msgrid,gps,strainmap.lonlims,strainmap.latlims,clim_straim,borders,places,...
        'Max shear (/yr)', faults, 0)%, poly)
    % colormap(ax, acton)
    colormap(ax, lajolla)
    
    ax = nexttile; hold on
    % plt_data(lon,lat,i2grid,gps,strainmap.lonlims,strainmap.latlims,clim_straim,borders,places,...
    %     'Second invariant of strain rate (/yr)', faults)%, poly)
    %JF
    plt_data(lon_grid,lat_grid,i2grid,gps,strainmap.lonlims,strainmap.latlims,clim_straim,borders,places,...
        'Second invariant of strain rate (/yr)', faults, 0)%, poly)
    % colormap(ax, acton)
    colormap(ax, lajolla)
    
    ax = nexttile; hold on
    % plt_data(lon,lat,digrid,gps,strainmap.lonlims,strainmap.latlims,[-1e-7, 1e-7],borders,places,...
    %     'Dilatation (/yr)', faults)%, poly)
    %JF
    plt_data(lon_grid,lat_grid,digrid,gps,strainmap.lonlims,strainmap.latlims,[-1e-7, 1e-7],borders,places,...
        'Dilatation (/yr)', faults, 0)%, poly)
    colormap(ax, bam)
    
    %% plot strains
    
    % nexttile; hold on
    % plt_data(lon,lat,ssgrid,gps,strainmap.lonlims,strainmap.latlims,clim,borders,places,...
    %     'Savage and Simpson [1997] function of strain rate (/yr)')
    % colormap(acton)
    
    sgtitle(strrep(outdir, '_', ' ')) 
    
    %saveas(gcf,outdir+'/vel_strain.png');
    saveas(gcf,strcat(outdir,'vel_strain.png'));
end 
%% if second input exists, plot differences
if plot_diff
    %%%%%
    % East diff
    egrid_diff = egrid - egrid2;
    data_flat_ediff = egrid_diff(:);
    data_flat_ediff = data_flat_ediff(~isnan(data_flat_ediff));
    clim_ediff = round([prctile(data_flat_ediff, 1) prctile(data_flat_ediff, 99)]);

    % North diff
    ngrid_diff = ngrid - ngrid2;
    data_flat_ndiff = ngrid_diff(:);
    data_flat_ndiff = data_flat_ndiff(~isnan(data_flat_ndiff));
    clim_ndiff = round([prctile(data_flat_ndiff, 1) prctile(data_flat_ndiff, 99)]);
    clim_ndiff = [-1,1];
    % Up diff
    ugrid_diff = ugrid - ugrid2;
    data_flat_udiff = ugrid_diff(:);
    data_flat_udiff = data_flat_udiff(~isnan(data_flat_udiff));
    clim_udiff = round([prctile(data_flat_udiff, 1) prctile(data_flat_udiff, 99)]);
    
    % Strain differences
    msgrid_diff = msgrid - msgrid2;
    i2grid_diff = i2grid - i2grid2;
    digrid_diff = digrid - digrid2;

    %%%%%
    
    % create new figure for differences
    figure()
    tiledlayout(2,3,'TileSpacing','compact')
    axis equal;
    
    % East velocity difference
    ax = nexttile; hold on
    plt_data(lon_grid, lat_grid, egrid_diff, gps, strainmap.lonlims, strainmap.latlims, clim_ndiff, borders, places, 'East vel diff (mm/yr)', faults, 0)
    colormap(ax, vik)

    % North velocity difference
    ax = nexttile; hold on
    plt_data(lon_grid, lat_grid, ngrid_diff, gps, strainmap.lonlims, strainmap.latlims, clim_ndiff, borders, places, 'North vel diff (mm/yr)', faults, 0)
    colormap(ax, vik)

    % Up velocity difference
    ax = nexttile; hold on
    plt_data(lon_grid, lat_grid, ugrid_diff, gps, strainmap.lonlims, strainmap.latlims, clim_ndiff, borders, places, 'Up vel diff (mm/yr)', faults, 0)
    colormap(ax, vik)

    % Max shear strain diff
    ax = nexttile; hold on
    plt_data(lon_grid, lat_grid, msgrid_diff, gps, strainmap.lonlims, strainmap.latlims, clim_straim, borders, places, 'Max shear diff (/yr)', faults, 0)
    colormap(ax, lajolla)

    % Second invariant diff
    ax = nexttile; hold on
    plt_data(lon_grid, lat_grid, i2grid_diff, gps, strainmap.lonlims, strainmap.latlims, clim_straim, borders, places, 'Second inv diff (/yr)', faults, 0)
    colormap(ax, lajolla)

    % Dilatation diff
    ax = nexttile; hold on
    plt_data(lon_grid, lat_grid, digrid_diff, gps, strainmap.lonlims, strainmap.latlims, [-1e-7 1e-7], borders, places, 'Dilatation diff (/yr)', faults, 0)
    colormap(ax, bam)

    sgtitle(['Difference: ' strrep(outdir, '_', ' ') ' - ' strrep(outdir2, '_', ' ')])
    
    % save figure
    saveas(gcf, strcat(outdir, 'vel_strain_diff.png'));
end


    % ,asc_poly, dsc_poly, faultfile
    %% nested plotting functions -----------------------------------------------------
    function plt_data(lon, lat, data, gps, lonlim, latlim, clim, borders, places, titlestr, faults, gps_plot)
    
    % figure(); hold on
    
    % plot input data
    imagesc(lon,lat,data);
    
    % % add country borders
    % for ii = 1:length(places)
    %     b_ind = find(strcmp(borders.places,places(ii)));
    %     plot(borders.lon{b_ind},borders.lat{b_ind},'k')
    % end
    
    % plot faults
    plot(faults(:,1),faults(:,2),'k')
    
    % plot insar coverage
    %plot(poly(:,1),poly(:,2),'w')
    
    
    % plot gps arrows
    if gps_plot && ~isempty(gps)
        try
            for igps = 1:length(gps)
                gpssite = gps(igps).site;
                if isempty(gpssite)
                    continue
                end
                lon = [gpssite.lon]';
                lat = [gpssite.lat]';
                east = cellfun(@(v) v(1), {gpssite.vel})';
                north = cellfun(@(v) v(2), {gpssite.vel})';
                sc = 0.025;
                clat = mean(latlim);  % auto center latitude based on map
                quiver(lon, lat, sc/cosd(clat)*east, sc*north, 0, 'color', [0.5 0.5 0.5], 'linewidth', 0.5);
           
            end
        catch ME
            warning('Problem plotting GPS data: %s', ME.message);
        end
    end
    
    axis xy
    xlim(lonlim)
    ylim(latlim)
    
    colorbar
    caxis(clim)
    
    
    title(titlestr)
    
    end
end