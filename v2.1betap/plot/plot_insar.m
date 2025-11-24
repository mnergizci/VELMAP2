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

function plot_insar(insarfit,gps,outdir, pngfile)

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

% pre-allocate
x = cell(1,length(insarfit));
y = cell(1,length(insarfit));
datamap = cell(size(x));
resimap = cell(size(x));
refmap = cell(size(x));
stackmap = cell(size(x));
ratemap = cell(size(x));
orbmap = cell(size(x));
atmmap = cell(size(x));
insarboundary = cell(size(x));
passdir = cell(size(x));

for ii = 1:length(insarfit)
    
    % lon and lat coord vectors
    x{ii} = insarfit(ii).ifghdr.xfirst ...
        + [0:(insarfit(ii).ifghdr.width-1)].*insarfit(ii).ifghdr.xstep;
    y{ii} = insarfit(ii).ifghdr.yfirst ...
        + [0:(insarfit(ii).ifghdr.length-1)].*insarfit(ii).ifghdr.ystep;
    
    % boundary of each frame
    [xx,yy] = meshgrid(x{ii},y{ii});
    xx = xx(:); yy = yy(:);
    xx(isnan(insarfit(ii).ratemap(:))) = []; yy(isnan(insarfit(ii).ratemap(:))) = [];
    
    boundaryind = boundary(xx(:),yy(:));
    insarboundary{ii} = [xx(boundaryind) yy(boundaryind)];
    
    passdir{ii} = insarfit(ii).ifghdr.passdir;
    
    datamap{ii} = insarfit(ii).stackmap + insarfit(ii).resmap;
%     stackresorb{ii} = stackres{ii} - insarfit(ii).orbmap;
    refmap{ii} = insarfit(ii).ratemap + insarfit(ii).resmap;
    resimap{ii} = insarfit(ii).resmap;
    stackmap{ii} = insarfit(ii).stackmap;
    ratemap{ii} = insarfit(ii).ratemap;
    orbmap{ii} = insarfit(ii).orbmap;
    atmmap{ii} = insarfit(ii).stackmap - insarfit(ii).orbmap - insarfit(ii).ratemap;
end

%% get asc and desc

asc_ind = strcmp(passdir,'A');
desc_ind = strcmp(passdir,'D');

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

% if pngfile == "model_los.png"
%     clim_asc = [-100 100];
%     clim_desc = [-100 100];
% end 

% clim_asc = [-20 20];
% clim_desc = [-20 20];

%% axes limits

lonlim = [min(cell2mat(x)) max(cell2mat(x))];
latlim = [min(cell2mat(y)) max(cell2mat(y))];

%% plot all together

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.2, 0.6, 0.6]);

tiledlayout(4,3,'TileSpacing','compact')
axis equal;

%%%%ascending section
nexttile; hold on
plot_data(x(asc_ind),y(asc_ind),refmap(asc_ind),insarboundary(asc_ind),...
    gps,lonlim,latlim,clim_asc,borders,places,'Referenced Asc')

nexttile; hold on
plot_data(x(asc_ind),y(asc_ind),datamap(asc_ind),insarboundary(asc_ind),...
    gps,lonlim,latlim,clim_asc,borders,places,'Asc')

% c=colorbar;
% c.Label.String = 'asc mm/yr';

nexttile; hold on
plot_data(x(desc_ind),y(desc_ind),refmap(desc_ind),insarboundary(desc_ind),...
    gps,lonlim,latlim,clim_desc,borders,places,'Referenced Dsc')


nexttile; hold on
plot_data(x(desc_ind),y(desc_ind),datamap(desc_ind),insarboundary(desc_ind),...
    gps,lonlim,latlim,clim_desc,borders,places,'Dsc')
% c=colorbar;
% c.Label.String = 'dsc mm/yr';


nexttile; hold on
plot_data(x(asc_ind),y(asc_ind),ratemap(asc_ind),insarboundary(asc_ind),...
    gps,lonlim,latlim,clim_asc,borders,places,'Model Asc')

nexttile; hold on
plot_data(x(asc_ind),y(asc_ind),orbmap(asc_ind),insarboundary(asc_ind),...
    gps,lonlim,latlim,clim_asc,borders,places,'Model Ramps')
c=colorbar;
c.Label.String = 'asc mm/yr';

nexttile; hold on
plot_data(x(desc_ind),y(desc_ind),ratemap(desc_ind),insarboundary(desc_ind),...
    gps,lonlim,latlim,clim_desc,borders,places,'Model Dsc')

nexttile; hold on
plot_data(x(desc_ind),y(desc_ind),orbmap(desc_ind),insarboundary(desc_ind),...
    gps,lonlim,latlim,clim_desc,borders,places,'Model Ramps')
c=colorbar;
c.Label.String = 'dsc mm/yr';

nexttile; hold on
plot_data(x(asc_ind),y(asc_ind),resimap(asc_ind),insarboundary(asc_ind),...
    gps,lonlim,latlim,clim_asc,borders,places,'Residual Asc')

nexttile; hold on
plot_data(x(asc_ind),y(asc_ind),atmmap(asc_ind),insarboundary(asc_ind),...
    gps,lonlim,latlim,clim_asc,borders,places,'Model Atm')
% c=colorbar;
% c.Label.String = 'asc mm/yr';

nexttile; hold on
plot_data(x(desc_ind),y(desc_ind),resimap(desc_ind),insarboundary(desc_ind),...
    gps,lonlim,latlim,clim_desc,borders,places,'Residual Dsc')

nexttile; hold on
plot_data(x(desc_ind),y(desc_ind),atmmap(desc_ind),insarboundary(desc_ind),...
    gps,lonlim,latlim,clim_desc,borders,places,'Model Atm')
% c=colorbar;
% c.Label.String = 'dsc mm/yr';
colormap(vik)

sgtitle(strrep(outdir, '_', ' ')) 

%saveas(gcf, outdir+'/los.png');
saveas(gcf, fullfile(outdir, pngfile));



%% nested plotting function
    function plot_data(x, y, data, insarboundary, gps, lonlim, latlim, clim, borders, places, titlestr)
        hold on
        for jj = 1:length(x)
            imagesc(x{jj}, y{jj}, data{jj}, 'AlphaData', ~isnan(data{jj}))
        end

        % Country borders (optional)
        for ii = 1:length(places)
            b_ind = find(strcmp(borders.places, places(ii)));
            % plot(borders.lon{b_ind}, borders.lat{b_ind}, 'r')
        end

        caxis(clim)
        xlim(lonlim)
        ylim(latlim)
        title(titlestr)
    end

end