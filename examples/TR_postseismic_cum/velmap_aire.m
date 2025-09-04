function velmap_aire(cfgfile)
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
clearvars -except cfgfile
clc
tic

%% MN setup will be read in velmap('config.file'). 
% %% 0. setup
% % input config file
% cfgfile = 'velmap_cum.conf';

%% Define VELMAP root path
VELMAP_path = '/scratch/eemne/github/VELMAP/';
% Add paths to required VELMAP subdirectories
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'pilib')));
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'static')));
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'mesh')));
addpath(genpath(fullfile(VELMAP_path, 'v2.1betap', 'plot')));

fprintf('readparfile -> %s\n', which('readparfile'));
fprintf('getparval   -> %s\n', which('getparval'));

%% validate & show cfgfile
if nargin < 1 || isempty(cfgfile)
    error('velmap_aire:NoConfig','cfgfile input is required');
end
if ~ischar(cfgfile) && ~isstring(cfgfile)
    error('velmap_aire:BadType','cfgfile must be char or string');
end
cfgfile = char(cfgfile);  % ensure char for older code
fprintf('cfgfile = %s\n', cfgfile);
fprintf('pwd     = %s\n', pwd);
if ~isfile(cfgfile)
    error('velmap_aire:MissingFile','Config file not found: %s', cfgfile);
end


%% 1. read config file

fprintf('===> reading config file ...\n');
% 1 = legacy mode, 0 = new parameter file setup
[par,gpspar,insarpar,smpar,tssmpar,profpar, sboipar] = readparfile(cfgfile,1, 1); %make sure (config, 1, 1) for SBOI other case should be (config, 1 ,0)

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

%% 3. prepare gps observations

% Save
[~, fname, ~] = fileparts(gpspar.filename{1});  % strip path + extension
gpsfile = [fname '.mat'];                       % new .mat filename

fprintf('===> loading gps data ...\n');
if par.reload_gps
    if gpspar.ngpsfile>0    
        [gps]=loadgps(gpspar);        %load gps data
        for i=1:gpspar.ngpsfile
            [gps(i).site]=tidygps(trim,gps(i).site); %remove gps outside of the mesh
            gps(i).nsite=length(gps(i).site);
        end
    else
        fprintf('===> no gps data to load...exiting...\n');
        return;
    end

else
    fprintf('===> load %s...\n', gpsfile);
    load(gpsfile, 'gps');
end


%% 4. prepare insar observations

if par.reload_insar
    if insarpar.ninsarfile>0
        if  insarpar.activate > 0 
            fprintf('===> loading insar data ...\n');
            [insar]=loadlics(insarpar);  %load insar data in lics format
        else
            fprintf('===>insar not activated to load...\n');
            insar=[];  
        end
    else
        fprintf('===> no insar data to load...\n');
        insar=[];

    end
    
else
    fprintf('===> load insar.mat');
    load insar.mat   
end 

% plottrim(trim,gps,insar);
% drawnow;
% saveas(gcf,'mesh.png');

%% 4.1 prepare sboi observations

if par.reload_sboi
    if  sboipar.activate > 0
        fprintf('===> loading sboi data ...\n');
        [sboi]=loadlics_sboi(sboipar);  %load insar data in lics format
        nsboi = length(sboi);
%         for i = 1:nsboi
%               sboi(i).azi = sboi(i).azi + 180;
%         end

    else
        fprintf('===> no sboi data to load...\n');
        sboi=[];

    end
    
else
    fprintf('===> load sboi.mat');
    load sboi.mat   
end 

%plotting data and mesh in together
% if par.reload_sboi
%     if sboipar.activate>0
%         plottrim(trim,gps,insar,sboi);
%         drawnow;
%         saveas(gcf,'mesh_with_boi.png');
%     end
% end

%% 5. find best smoothing factor

if smpar.smf==0
  fprintf('===>  processing for all smoothing factors...\n');
  smpar.smf = (smpar.smf_min:smpar.smf_int:smpar.smf_max);
  smpar.smf = 10.^smpar.smf;
  
elseif smpar.smf==999 % currently non-functioning
  fprintf('===>  find the best smoothing factor...\n');
  smpar.smf = lcurve_vmp(trim,smpar,gps,insar,1,par.outdir);
  
end


%% solution for each smoothing factor

%update invenu
invenu=getinvenu(gps,insar);
nsmf=length(smpar.smf);

for i=1:nsmf
  fprintf('===> processing smoothing factor %d/%d\n',i,nsmf);
  %output directory
  orb_str = sprintf('orb%d-%d', insarpar.orbdegree, sboipar.orbdegree);
  atm_str = sprintf('atm%d-%d', insarpar.atmdegree, sboipar.atmdegree);
  lk_str  = sprintf('lk%d', insarpar.lksx);
  mesh_str = sprintf('mesh%.2f', par.mesh_dx);
  smf_str  = sprintf('smf%+4.2f', log10(smpar.smf(i)));
  n_sboi = length(sboi);
  n_insar = length(insar);
  n_gps = gps(1).nsite;% + gps(2).nsite;

%   outdir = sprintf('%s_%s_%s_insars%d_sbois%d_gps%d_%s_%s_%s_%s/', ...
%     char(par.outdir), smf_str,n_insar,n_sboi, n_gps, ...
%     orb_str, atm_str, lk_str, mesh_str);
  outdir = sprintf('%s_%s_insars%d_sbois%d_gps%d_%s_%s_%s_%s/', ...
    char(par.outdir), smf_str, n_insar, n_sboi, n_gps, ...
    orb_str, atm_str, lk_str, mesh_str);

  if ~exist(outdir,'dir')
    mkdir(outdir)
  end
  save(fullfile(outdir, 'sboi.mat'), 'sboi')
  save(fullfile(outdir, 'gps.mat'), 'gps')
  save(fullfile(outdir, 'insar.mat'), 'insar')
  
  %%%%extras
  plottrim(trim,gps,insar);
  drawnow;
  saveas(gcf, fullfile(outdir, 'mesh.png'));
  
  if par.reload_sboi && sboipar.activate > 0 && exist('sboi','var')
    plottrim(trim, gps, insar, sboi);
    drawnow;
    % saveas(gcf, fullfile(outdir, 'mesh_with_boi.png'));
    saveas(gcf, fullfile(outdir, 'mesh_with_boi.png'));
  end
  %%%%%%
  
  %% 6. solve system of equations

  %   [fitmodel,vcmmodel,wrss,rough]=solve_vmp_single(trim,smpar.smf(i),gps,insar,outdir,1); % single core version
  [fitmodel,vcmmodel,wrss,rough]=solve_vmp(trim,smpar.smf(i),gps,insar,sboi,outdir,1); % parallel version
  lcv=[smpar.smf(i) rough wrss log10(smpar.smf(i))]; %JRW add
  
  % JRW add - Andrew edited to avoid changing directory and change deprecated dlmwrite
  % to writematrix
  writematrix(lcv,[outdir 'lcv.dat'],'delimiter','\t')
  save([outdir 'fitmodel'],'fitmodel','-v7.3');
  save([outdir 'vcmmodel'],'vcmmodel','-v7.3');

  %% 7. forward calculation

  % output fitted velocity field
  fprintf('===> output fitted velocity field... \n');
  fitvtx = fitmodel2vel(trim,fitmodel,vcmmodel,invenu,outdir);
  
  save([outdir 'fitvtx'],'fitvtx','-v7.3');

  % forward calculation for gps
  fprintf('===> forward calculating... \n');
  gpsfit = gpsfwd(trim,fitmodel,vcmmodel,invenu,gps,outdir);
  
  save([outdir 'gpsfit'],'gpsfit','-v7.3');
  
  % forward calculation for insar
  if insarpar.ninsarfile>0
    if insarpar.activate > 0
        insarfit = insarfwd(insar,trim,fitmodel,invenu,outdir,gps);
        save([outdir 'insarfit'],'insarfit','-v7.3');    
    end
  end
  
  if insarpar.ninsarfile>0
      if sboipar.activate > 0
        sboifit = insarfwd(sboi,trim,fitmodel,invenu,outdir,gps, 1);
        save([outdir 'sboifit'],'sboifit','-v7.3');    
      end
  end
  %Save all what we create so far
  save([outdir 'precrash.mat'], '-v7.3');

  %% 8. strain rate

  % calculate strain rate
  fprintf('===> calculating strain rate... \n');
  nvtx=length(trim.x);
  fitvel=(reshape([fitvtx.vel],[],nvtx))';
%   [strain]=vel2strain_tri(trim,fitvel,outdir);
  [strain,eulervec]=vel2strain_savage(trim,fitvel,vcmmodel,1,outdir,2);
%   [strain,eulervec]=vel2strain_savage(trim,fitvel,vcmmodel,2,outdir,2);

  %% 9. profile
  %make profile for the velocity field
  if profpar.profflag==1
    fprintf('===> making profiles... \n');
    %interactively extract profile
    profdef=char('profdef.dat');
    if ~exist(profdef,'file')
      plotvel(fitvtx,gps,gpsfit,faults,outdir);
      extractprof(prof.swath,prof.step);
    end

    %read profile
    [prof]=profdefine(profdef);

    %extract fault position on the profile
    if ~exist('faultonprof.dat','file')
      extractfaultonprof(prof,faults);
    end

    %calculate profile
    nprof=length(prof);
    profdir=strcat(outdir,'prof/');
    if ~exist(profdir,'dir')
      mkdir(profdir);
    end
    for iprof=1:nprof
      fprintf('making profile %d/%d\n',iprof,nprof);
      if prof(iprof).swath==0
        make_profline_vel(trim,fitmodel,vcmmodel,invenu,prof(iprof),profdir);
      else
        make_profswath_vel(trim.x,trim.y,fitvel,vcmmodel(1:2*nvtx,1:2*nvtx),prof(iprof),profdir);
      end
      %make profile for the observed GPS data
      %low efficiency to make profile for each site once a time
      gpsprofdef=prof(iprof);
      gpsprof=[];
      gpsprofdef.swath=gpsswath;
      for igf=1:gpspar.ngpsfile
        ns=length(gps(i).site);
        for is=1:ns
          isite=gps(igf).site(is);
          igpsprof=make_profswath_vel(isite.lon,isite.lat,isite.vel,isite.vcm,gpsprofdef);
          gpsprof=[gpsprof;igpsprof];
        end
      end
      if size(gpsprof,1)>0
        outfile=strcat(profdir,gpsprofdef.id,'.prof_gps');
        save(char(outfile),'gpsprof','-ASCII');
      end

      %make profile for the observed InSAR data
      insarprofdef=prof(iprof);
      insarprofdef.swath=insarswath;
      for isf=1:insarpar.ninsarfile
        %stackmap=insar(isf).stackmap; 
        %using original resolution stackmap
        %ifghdr=rsc2hdr(char(strcat(insarpar.dir(isf),'ratemap/ifg.rsc')));
        %ifghdr=rsc2hdr(char(strcat(insarpar.dir(isf),'ifghdr.mat')));
        ifghdr=char(strcat(insarpar.dir(isf),'ifghdr.mat')); %old version of pi-rate 2.0 or earlier - JRW add
        load(ifghdr); %JRW add
        %stackmap=readmat(char(strcat(insarpar.dir(isf),'ratemap/stackmap.dat')),ifghdr.length,ifghdr.width,1);
        stackmap=readmat(char(strcat(insarpar.dir(isf),'stackmap.dat')),ifghdr.length,ifghdr.width,1); %JRW add
        [sarprof_pt,sarprof] = make_prof(stackmap,insarprofdef,ifghdr);
        if size(sarprof,1)>0
          sarprof_pt=double(sarprof_pt);
          sarprof=double(sarprof);
          outfile=strcat(profdir,num2str(insarprofdef.id),'.prof_insar',num2str(isf,'%02d'));
          save(char(outfile),'sarprof','-ASCII');
          outfile=strcat(profdir,num2str(insarprofdef.id),'.prof_insar',num2str(isf,'%02d'),'_pt');
          save(char(outfile),'sarprof_pt','-ASCII');
        end
      end
    end
  end

  %% 10. make velocity field on a regular grid
%   disp('===> NOT making velocity field on a regular grid...')
  %make_grid_vel(trim,fitmodel,vcmmodel,invenu,grdvel.dx,grdvel.dy,outdir);
  
  %% 11. plot results
  
    % dump insarfit
    if insarpar.ninsarfile~=0
        if insarpar.activate > 0
            plot_insar(insarfit,gps, outdir, 'model_los.png');
            dump_insarfit(insarfit, insar, outdir);
        end
        if sboipar.activate > 0
            plot_insar(sboifit, gps, outdir, 'model_azi.png');
            dump_insarfit(sboifit, sboi, outdir, 'sboifit');
        end
    end
    %plot_vel_strain
    myplot_vel_strain(outdir)
    myplot_vel_std(outdir)

    %% 12. Return full resolution LiCSBAS output
    %post_process_full_res(outdir)
end

fprintf('====finished successfully, congratulations!====\n');

toc
