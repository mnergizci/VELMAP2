function [] = plottrim(trim,gps,insar, sboi)
%===========================================
%function [] = plottrim(trim,gps,insar)
%
% Plot triangular mesh and insar data together
%
% Input:
%  trim: triangular mesh
%  gps:  gps data
%  insar: insar data
%  sboi: sboi data
%
% Output:
%  None
%
% Hua Wang @ Uni Leeds, 01/10/2009
%
% 23/09/2015 HW: revise for new gps structure
%                renamed as plotrim as plotmesh is a embeded function
%===========================================

%plot triangular mesh
figure

if nargin>2
  ninsar=length(insar);
  for is=1:ninsar
    [rows,cols]=size(insar(is).stackmap);
    n=rows*cols;
    vrate=reshape(insar(is).stackmap',n,1);
    [xx,yy]=meshgrid(1:cols,1:rows);
    xxv=reshape(xx',n,1);
    yyv=reshape(yy',n,1);
    xxv = insar(is).ifghdr.xfirst+(xxv-1)*insar(is).ifghdr.xstep;
    yyv = insar(is).ifghdr.yfirst+(yyv-1)*insar(is).ifghdr.ystep;
    clear('xx','yy');
    xxv(isnan(vrate))=[];
    yyv(isnan(vrate))=[];
    vrate(isnan(vrate))=[];
    r=is/ninsar;
    g=1-is/ninsar;
    b=1-is/ninsar;
    plot(xxv,yyv,'.','MarkerSize',10,'MarkerEdgeColor',[r g b]);
    hold on
  end
end

if nargin>3
  nsboi = length(sboi);
  for is = 1:nsboi
    [rows, cols] = size(sboi(is).stackmap);
    n = rows * cols;
    vrate = reshape(sboi(is).stackmap', n, 1);
    [xx, yy] = meshgrid(1:cols, 1:rows);
    xxv = reshape(xx', n, 1);
    yyv = reshape(yy', n, 1);
    xxv = sboi(is).ifghdr.xfirst + (xxv - 1) * sboi(is).ifghdr.xstep;
    yyv = sboi(is).ifghdr.yfirst + (yyv - 1) * sboi(is).ifghdr.ystep;
    xxv(isnan(vrate)) = [];
    yyv(isnan(vrate)) = [];
    vrate(isnan(vrate)) = [];
    % r = 0.2; g = 0.2; b = is/nsboi; % Light blue for sboi
    colors = hsv(nsboi);   % Vivid and colorful
    color = colors(mod(is-1, size(colors,1)) + 1, :);
    plot(xxv, yyv, '.', 'MarkerSize', 8, 'MarkerEdgeColor', color);
    hold on
  end
end  


triplot(trim.tri,trim.x,trim.y);

hold on


if nargin>1
  scale=50;
  for igf = 1:length(gps)
    site = gps(igf).site;
    ngps = length(site);
    vel = zeros(ngps, 2);
    for i = 1:ngps
      if isfield(site(i), 'vel')
        vel(i,:) = site(i).vel(1:2);
      elseif isfield(site(i), 'tsvel')
        vel(i,:) = site(i).tsvel(1,1:2);
      end
    end
    quiver([site.lon]', [site.lat]', vel(:,1)/scale, vel(:,2)/scale, 'k', 'AutoScale', 'off');
  end
end


axis equal
axis([min(trim.x)-1,max(trim.x)+1,min(trim.y)-1,max(trim.y)+1]);