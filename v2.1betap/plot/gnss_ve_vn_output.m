%=================================================================
% gnss_ve_vn_output.m
% Plot velocities and strains from velmap output.
% Tiledlayout requires 2019b or later.
%
% Muhammet Nergizci @ leeds, 26/04/2025
%=================================================================

function gnss_ve_vn_output(outdir, gps)

% Handle optional inputs
if nargin < 2 || isempty(gps)
    if exist(fullfile(outdir, 'gps.mat'), 'file') == 2
        S = load(fullfile(outdir, 'gps.mat'));
        gps = S.gps;
    else
        warning('gps.mat not found. GPS data will not be plotted.');
        gps = [];
    end
end


% Define input file paths
velfile     = fullfile(outdir, 'velfit.dat');
strainfile  = fullfile(outdir, 'strain_savage_nring1.dat');
faultfile   = 'faults/gem_active_faults.gmt';

% Read velocity and strain data
vel    = readmatrix(velfile);
strain = readmatrix(strainfile);

% === Step 2: Extract lon, lat, E, N ===
x = vel(:, 1);   % x
y = vel(:, 2);   % y
E   = vel(:, 3);   % east velocity
N   = vel(:, 4);   % north velocity
sE = vel(:,6);
sN = vel(:,7);

% Filter NaNs
valid = all(~isnan([x y E N sE sN]), 2);
x  = x(valid);
y  = y(valid);
E  = E(valid);
N  = N(valid);
sE = sE(valid);
sN = sN(valid);

% Store into struct
gnss_field = struct('x', x, 'y', y, 'E', E, 'N', N, 'sE', sE, 'sN', sN);

%%saving
save('gnss_decomp.mat', 'gnss_field');

% You may now use gnss.mat in velmap tools or plotting scripts
disp('GNSS struct saved to gnss.mat');



end