%addpath(genpath('v2.1betap'))

%load precrash.mat;
[insar]=loadlic1lk(insarpar);
%save('insar_1lk','insar','-v7.3');
%smfdir = 'out/'
%load(strcat(smfdir, 'fitmodel.mat'));
%[insarfit2]=insarfwd(insar,trim,fitmodel,invenu,smfdir,gps);
[insarfit2]=insarfwd(insar,trim,fitmodel,invenu,outdir,gps);
%cd (smfdir)
%save('insarfit2','insarfit2','-v7.3');
%cd ../

outputgeotiffs

fprintf('====Finished successfully. Saving insarfit2.mat====\n');
cd (outdir)
save('insarfit2','insarfit2','-v7.3');
cd ../
