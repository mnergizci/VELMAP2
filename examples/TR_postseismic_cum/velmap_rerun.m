%=================================================================
% velmap
% Main script to make a crustal velocity map using geodetic data.
%
% Calls functions from the main, pilib, and plotting directories.
% Written for Matlab 2019a.
%
% Can use the original parameter file setup or a new version. Toggle for
% this is found in readparfile.
%
% Original - Hua Wang @ Leeds, 17/09/2009
% Updated - Andrew Watson @ leeds, 15/06/2021
% Updated - Jin Fang @ Leeds, 7/3/2025
%=================================================================
close all
clear;clc
tic
%% 0. setup
% input config file
cfgfile = 'velmap.conf';

% Define VELMAP root path
VELMAP_path = 'C:\Users\mmuha\OneDrive - University of Leeds\1.Phd_project\1.velmap_arie\VELMAP';
VELMAP_path = '/mnt/scratch/eemne/VELMAP/'
VELMAP_path = '/nfs/a1/eemne/VELMAP'
% Add paths to required VELMAP subdirectories
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'pilib')));
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'static')));
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'mesh')));
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'plot')));
outdir='output_vel_mskd-withnan0.2_smf-1.20_no_adding180_insars0_sbois24_orb2-1_atm1-0_lk5_mesh0.20/'

%% 1. read config file

fprintf('===> reading config file ...\n');
% 1 = legacy mode, 0 = new parameter file setup
[par,gpspar,insarpar,smpar,tssmpar,profpar, sboipar] = readparfile(cfgfile,1, 1); %make sure (config, 1, 1) for SBOI other case should be (config, 1 ,0)
load([outdir 'gps.mat']);
load([outdir 'gpsfit.mat']);
load([outdir 'insar.mat']);
load([outdir 'insarfit.mat']);
load([outdir 'fitmodel.mat']);
load([outdir 'fitvtx.mat']);
load([outdir 'sboi.mat']);
load([outdir 'sboifit.mat']);
load([outdir 'vcmmodel.mat']);


%% 2. make triangular mesh
% fault location, range of study area needed

fprintf('===> reading mesh ...\n');

if exist(par.meshfile,'file')
    if strcmp(par.meshfile(end-3:end),'.msh')
        trim=gid2mat(par.meshfile);
        
    elseif strcmp(par.meshfile(end-3:end),'.mat')
        %makemesh(par.meshfile, par.mesh_dx, par.mesh_dy)
        load(par.meshfile) 
    end
    
else
    fprintf('mesh file is not available, lets create a new one');
    makemesh('TR_EAF_mesh.mat', 0.2, 0.2); %%play with inside
    load('TR_EAF_mesh.mat');

end
  %% 8. strain rate

  % calculate strain rate
  fprintf('===> calculating strain rate... \n');
  nvtx=length(trim.x);
  fitvel=(reshape([fitvtx.vel],[],nvtx))';
%   [strain]=vel2strain_tri(trim,fitvel,outdir);
  [strain,eulervec]=vel2strain_savage(trim,fitvel,vcmmodel,1,outdir,2);
%   [strain,eulervec]=vel2strain_savage(trim,fitvel,vcmmodel,2,outdir,2);

  %% 10. make velocity field on a regular grid
%   disp('===> NOT making velocity field on a regular grid...')
  %make_grid_vel(trim,fitmodel,vcmmodel,invenu,grdvel.dx,grdvel.dy,outdir);
  
  %% 11. plot results
  
    % dump insarfit
    if insarpar.ninsarfile~=0
        if insarpar.activate > 0
%             plot_insar(insarfit,gps, outdir, 'model_los.png');
%             plot_insar_double(insarfit,gps, outdir, 'model_los.png');
%             plot_insar_fourpanel(insarfit,gps, outdir, 'model_los.png');
%             dump_insarfit(insarfit, insar, outdir);
        end
        if sboipar.activate > 0
            plot_insar(sboifit, gps, outdir, 'model_azi.png');
            plot_insar_double(sboifit, gps, outdir, 'model_azi.png');
            plot_insar_fourpanel(sboifit, gps, outdir, 'model_azi.png');
            dump_insarfit(sboifit, sboi, outdir, 'sboifit');
        end
    end
    %plot_vel_strain
    myplot_vel_strain(outdir)
    myplot_vel_std(outdir)

    %% 12. Return full resolution LiCSBAS output
    %post_process_full_res(outdir)
fprintf('====finished successfully, congratulations!====\n');

toc
