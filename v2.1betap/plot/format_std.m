function [egrid,ngrid,ugrid] = format_std(vel,meshspacing)
%=================================================================
% format_strain.m
% Format vel outputs into plottable grids.
%
% Andrew Watson @ leeds, 17/02/2022
% Muhammet Nergizci @ leeds
%=================================================================

% gen uniform mesh
lon = vel(:,1); lat = vel(:,2);
[longrid,latgrid] = meshgrid(min(lon):meshspacing:max(lon),...
    min(lat):meshspacing:max(lat));

% interpolate onto mesh
einterp = scatteredInterpolant(lon,lat,vel(:,6),'linear','none');
egrid = einterp(longrid,latgrid);

ninterp = scatteredInterpolant(lon,lat,vel(:,7),'linear','none');
ngrid = ninterp(longrid,latgrid);

uinterp = scatteredInterpolant(lon,lat,vel(:,8),'linear','none');
ugrid = uinterp(longrid,latgrid);

end

