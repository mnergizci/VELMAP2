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
load([outdir 'fitmodel.mat']);
load([outdir 'vcmmodel.mat']);
load([outdir 'gps.mat']);
load([outdir 'insar.mat']);


%% 2. make triangular mesh
% fault location, range of study area needed
fprintf('mesh file is not available, lets create a new one');
makemesh('TR_EAF_meshhigh.mat', 0.01, 0.01); %%play with inside
load('TR_EAF_meshhigh.mat');


%% 3. strain rate
invenu=getinvenu(gps,insar);
% output fitted velocity field
fprintf('===> output fitted velocity field... \n');
fitvtx = fitmodel2vel(trim,fitmodel,vcmmodel,invenu,outdir);

save([outdir 'fitvtx'],'fitvtx','-v7.3');

