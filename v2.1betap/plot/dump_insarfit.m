function dump_insarfit(insarfit,insar, outdir, prefix)
%=================================================================
% dump_insarfit.m
% Export insarfit structures into XYZ files for GMT plotting.
%
% Dr. Jin Fang, Dr. Qi Ou, Muhammet Nergizci @ Leeds, 27/04/2025
%=================================================================

if nargin < 3
    outdir = '.'; % default: current folder
end

if nargin < 4
    prefix = 'insarfit'; % default
end


ninsar = length(insarfit);

for i=1:ninsar
  disp(['Working on ',prefix,' ',insar(i).proc.ID]) 
  stackmap=insarfit(i).stackmap+insarfit(i).resmap;
  a(:,:,1)=stackmap;             %observed InSAR
  a(:,:,2)=insarfit(i).stackmap; % Modeled InSAR
  a(:,:,3)=insarfit(i).resmap;   % Residual
  a(:,:,4)=insarfit(i).ratemap;  % Interpolated rate from mesh without orbit and atm
  a(:,:,5)=insarfit(i).orbmap;   % Orbit correction
  
  % Grid size 
  [rows, cols] = size(stackmap);
  n = rows * cols;
  
  % Flatten arrays
  [xx, yy] = meshgrid(1:cols, 1:rows);
  xxv = reshape(xx', n, 1);
  yyv = reshape(yy', n, 1);
  xxv = insarfit(i).ifghdr.xfirst + (xxv-1) * insarfit(i).ifghdr.xstep;
  yyv = insarfit(i).ifghdr.yfirst + (yyv-1) * insarfit(i).ifghdr.ystep;

  a1=reshape(a(:,:,1)',n,1);
  a2=reshape(a(:,:,2)',n,1);
  a3=reshape(a(:,:,3)',n,1);
  a4=reshape(a(:,:,4)',n,1);
  a5=reshape(a(:,:,5)',n,1);

  % Remove NaNs based on stackmap
  mask = ~isnan(a1);
  xxv = xxv(mask);
  yyv = yyv(mask);
  a1 = a1(mask);
  a2 = a2(mask);
  a3 = a3(mask);
  a4 = a4(mask);
  a5 = a5(mask);

  out=[xxv yyv a1 a2 a3 a4 a5];

  
  outfile = fullfile(outdir, [prefix, insar(i).proc.ID, '.xyz']);
  dlmwrite(outfile, out, 'precision','%4.4f','delimiter','\t')
  clear stackmap a xx yy xxv yyv a1 a2 a3 a4 a5 out mask
  
end

end