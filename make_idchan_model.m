function [] = make_idchan_model(mhke_x,mhke_y,flood)
% generate an idealized channel setup 
%
% Dependencies - 
%   a.) fvcom preprocessing toolbox (http://code.google.com/p/fvcom-toolbox/)
%   b.) gmsh (http://geuz.org/gmsh/), note they have binaries and the OSX ones work fine
%   d.) matlab netcdf toolbox (see instructions here to acquire:  http://code.google.com/p/fvcom-toolbox/wiki/Installation)
%   e.) others - google first, if it doesn't exist on the web, contact Geoff.
%   
%
%  - Geoff Cowles
%
% set paths and problem parameters

optparams

%addpath(genpath('/opt/matlab/fvcom-toolbox'));
%addpath('/opt/matlab/Mesh2Dv24');  
%clear all; close all;

%path_to_gmsh = '/usr/local/bin/gmsh';

%--------------------------------------------------------------------------
% model - specific setups this is the user-defined section
%--------------------------------------------------------------------------
casename = 'ic';
modelid = '01';


%background_res = 100.;   
%island_res = 100;     
%open_boundary_res = 1000;    
%min_res = min([background_res,island_res,open_boundary_res]);
min_res = 20;
gmsh_args   = ['-v 0 -2 -algo frontal -clmin '  num2str(min_res)];


make_node_lists = 0;
vecres = 100; %resolution of structured grid at which to plot vectors

% set geometry parameters
%chan_length   = 5000;
%chan_width    = 1000;
%fillet_radius = 500;
%ob_radius     = 5000;
%ds_channel    = 150; %mesh spacing in the channel
%ds_outer      = 500.; %mesh spacing at the open boundary

% bathymetry
%chan_depth      = 10.; %depth in the channel
%max_depth       = 75.; %depth at open boundary
%bath_smooth_its = 5;

% sponge damping coefficient and radius
%SpongeCoeffInflow = 1e-4;
%SpongeCoeffOutflow = 1e-3;
%SpongeRadInflow   = 2*ds_outer;
%SpongeRadOutflow   = 2*ds_outer;


%MHKE
%mhke_x = [];%-2000,-2000,2000];% turbine x-position
%mhke_y = [];%-200,200,0];   % turbine y-position
nmhke = numel(mhke_x);
mhke_At = At*ones(nmhke,1); % turbine area
mhke_Cp = Cp*ones(nmhke,1); % turbine power coefficient
%refine_mhke=1;
%dxR = 1000;
%dyR = 100;

% flow direction [=1, flood (positive x) = 0, ebb (negative x)]
%flood = 0;

% head (free surface elevation difference between inflow and outflow in m
%head = 0.5;

%--------------------------------------------------------------------------
% output file names 
%--------------------------------------------------------------------------

obcfile = [casename modelid '_obc.dat']; delete(obcfile);
corfile = [casename modelid '_cor.dat']; delete(corfile);
depfile = [casename modelid '_dep.dat']; delete(depfile);
spgfile = [casename modelid '_spg.dat']; delete(spgfile);
grdfile = [casename modelid '_grd.dat']; delete(grdfile);
elvfile = [casename modelid '_elj.nc']; delete(elvfile);
mhkfile = [casename modelid '_mhke.dat']; delete(mhkfile);

%--------------------------------------------------------------------------
% set boundary conditions 
%--------------------------------------------------------------------------
if(flood)
open_boundary_bc_outflow  = 41; 
open_boundary_bc_inflow   = 42; 
else
  open_boundary_bc_outflow  = 42; 
open_boundary_bc_inflow   = 41; 
end;
  

%--------------------------------------------------------------------------
% make geo file for gmsh
%--------------------------------------------------------------------------
fname = 'ic.geo';
delete(fname);
system('cp ic.geo_template ic.geo');
find_and_replace(fname,'HALF_CHAN_WIDTH',num2str(0.5*chan_width));
find_and_replace(fname,'HALF_CONSTRICTION_WIDTH',num2str(0.5*constriction*chan_width));
find_and_replace(fname,'FILLET',num2str(fillet_radius));
find_and_replace(fname,'HALF_CHAN_LENGTH',num2str(0.5*chan_length));
find_and_replace(fname,'OB_RADIUS',num2str(ob_radius));
find_and_replace(fname,'DS_CHANNEL',num2str(ds_channel));
find_and_replace(fname,'DS_OUTER',num2str(ds_outer));

%--------------------------------------------------------------------
% initial mesh using gmsh and read back in
%--------------------------------------------------------------------
system([path_to_gmsh  ' ic.geo ' gmsh_args ' -o tmp.msh']);
fprintf('reading initial mesh back into matlab\n');
global Mobj;
Mobj = read_gmsh_mesh('msh','tmp.msh'); 

%--------------------------------------------------------------------
% refine MHKE zones 
%--------------------------------------------------------------------
if(refine_mhke==1);
  
  % save old x,y for calculating ob nodes
  xold = Mobj.x;
  yold = Mobj.y;
  
  % calculate xc/yc
  Mobj.xc = zeros(Mobj.nElems,1);
  Mobj.yc = zeros(Mobj.nElems,1);
  for i=1:Mobj.nElems
    nds = Mobj.tri(i,:);
    Mobj.xc(i) = sum(Mobj.x(nds))/3.;
    Mobj.yc(i) = sum(Mobj.y(nds))/3.;
  end;

  % mark all cells within mhke_boxes
  mark = false(Mobj.nElems,1);
  for m=1:numel(mhke_x);
    xb = [mhke_x(m)-dxR,mhke_x(m)+dxR,mhke_x(m)+dxR,mhke_x(m)-dxR];
    yb = [mhke_y(m)-dyR,mhke_y(m)-dyR,mhke_y(m)+dyR,mhke_y(m)+dyR];
    ins = inpolygon(Mobj.xc,Mobj.yc,xb,yb); mark(find(ins==1))=1;
  end;

  % refine using refine_model
  Mobj = refine_model(Mobj,mark,100,.02);
  
  % redo the open boundary node list, hopefully we do not introduce
  % split elements on the open boundary
  for o=1:Mobj.nObs
    nlist = Mobj.obc_nodes(o,1:Mobj.nObcNodes(o));
    for n=1:numel(nlist);
      xp = xold(nlist(n));
      yp = yold(nlist(n));
      [radi,ipt] = min(sqrt( (Mobj.x-xp).^2 + (Mobj.y-yp).^2));
      Mobj.obc_nodes(o,n) = ipt;
    end;
  end;
end;



[Mobj] = setup_metrics(Mobj);




%--------------------------------------------------------------------
% Make sure triangle connectivity is CCW 
%--------------------------------------------------------------------
nds = Mobj.tri(1,:);
xck = Mobj.x(nds);
yck = Mobj.y(nds);
area = (xck(2)-xck(1))*(yck(3)-yck(1)) - (xck(3)-xck(1))*(yck(2)-yck(1));
if(area < 0);
  ii = Mobj.tri(:,2);
  Mobj.tri(:,2) = Mobj.tri(:,3);
  Mobj.tri(:,3) = ii;
end;
clear ii;

%--------------------------------------------------------------------
% set latitude to zero for zero Coriolies 
%--------------------------------------------------------------------
Mobj.lon = zeros(Mobj.nVerts,1);
Mobj.lat = zeros(Mobj.nVerts,1);
Mobj = add_coriolis(Mobj); 

%--------------------------------------------------------------------
% interpolate bathymetry 
%--------------------------------------------------------------------
Mobj.h = zeros(Mobj.nVerts,1);
Mobj.have_bath = 1;
Mobj.h = chan_depth.*ones(Mobj.nVerts,1);
  
%--------------------------------------------------------------------
% set BC flags 
%--------------------------------------------------------------------
Mobj.obc_id = Mobj.obc_type;

%--------------------------------------------------------------------
% set analytical bathymetry 
%--------------------------------------------------------------------


Mobj.h = chan_depth*ones(Mobj.nVerts,1);
half_chan_length = 0.5*chan_length;
half_chan_width  = 0.5*chan_width;
for i=1:Mobj.nVerts
  xp = Mobj.x(i); 
  yp = Mobj.y(i);
  yloc = sign(yp)*min(abs(yp),half_chan_width);
  if(xp > half_chan_length)
    rad = sqrt( (xp-half_chan_length).^2 + (yp-yloc).^2);
    fac = rad/(.75*ob_radius);
    Mobj.h(i) = fac*max_depth + (1-fac)*chan_depth;
  elseif(xp < -half_chan_length)
    rad = sqrt( (xp+half_chan_length).^2 + (yp-yloc).^2);
    fac = rad/(.75*ob_radius);
    Mobj.h(i) = fac*max_depth + (1-fac)*chan_depth;
  end;
end;
Mobj.h = min(Mobj.h,max_depth);

% smooth it using explicit Poisson
Mobj.h = smoothfield(Mobj.h,Mobj,.5,bath_smooth_its);

%--------------------------------------------------------------------
% create a spatially-distributed bottom roughness lengthscale file  
% use a simple linear dependence on depth
%--------------------------------------------------------------------
%nc = netcdf(z0bfile,'clobber');
%nc.title = 'Example scheme for bottom roughness';
%nc('nele') = Mobj.nElems;
%nc{'z0b'} = ncfloat('nele');
%nc{'z0b'}.long_name = 'bottom roughness';
%nc{'z0b'}.units = 'm';
%hc = nodes2elems(Mobj.h,Mobj);  %bathy on the elements
%nc{'z0b'}(1:Mobj.nElems) =  .001*ones(Mobj.nElems,1); 
%ierr = close(nc);

%--------------------------------------------------------------------
% dump the Coriolis
%--------------------------------------------------------------------
write_FVCOM_cor(Mobj,corfile) ;

%--------------------------------------------------------------------
% dump mesh and connectivity 
%--------------------------------------------------------------------
write_FVCOM_grid(Mobj,grdfile);
write_FVCOM_bath(Mobj,depfile); 

%--------------------------------------------------------------------
% open boundary node list 
%--------------------------------------------------------------------
Mobj.obc_type(1:end) = 1; %tidal forcing
write_FVCOM_obc(Mobj,obcfile); 

%--------------------------------------------------------------------
% boundary elevation forcing 
%--------------------------------------------------------------------
MyTitle = 'idChan M2 test';
SpectralFile = elvfile;
Components = {'M2'}; % ,'S2','N2','K1','M4','O1'};
Period     = [44714.16]; %,43200.00,45570.05,86164.09,22357.08,92949.63];
Period     = [43200]; %exactly 12 hours for idealized sim
Period     = [-99];

nComps = numel(Period);
if(Mobj.nObs==0)
        warning('cannot setup spectral open boundary, there is no open boundary in the mesh struct')
        return
end;
Amp = zeros(sum(Mobj.nObcNodes),nComps);
Phase = zeros(sum(Mobj.nObcNodes),nComps);
cnt = 0;
for ob=1:Mobj.nObs
  nObcs = Mobj.nObcNodes(ob);
  for j=1:nObcs
    cnt = cnt + 1;
    ObcNodes(cnt) = Mobj.obc_nodes(ob,j);  %set the open boundary nodes
    for i=1:nObcs
      if(Mobj.obc_id(ob)==open_boundary_bc_inflow) % left side
        Amp(cnt,1:nComps) = head;
        Phase(cnt,1:nComps) = 0.0;
      else % right side
        Amp(cnt,1:nComps) = 0.0;
        Phase(cnt,1:nComps) = 0.0;
      end;
    end;
    
  end;
end;
EqAmp = zeros(nComps,1);
BetaL = EqAmp;
write_FVCOM_spectide(ObcNodes,Components,Period(1:nComps),Phase,Amp,EqAmp, ...
  BetaL,SpectralFile,MyTitle)

%--------------------------------------------------------------------
% add a sponge layer 
%--------------------------------------------------------------------
Mobj.nSponge = 0;
for ob=1:Mobj.nObs
  nObcs = Mobj.nObcNodes(ob);
  Mobj.nSponge = Mobj.nSponge + 1;
  Mobj.nSpongeNodes(Mobj.nSponge) = nObcs;
  Mobj.sponge_nodes(Mobj.nSponge,1:nObcs) = Mobj.obc_nodes(ob,1:nObcs);
  
  if(Mobj.obc_id(ob)==open_boundary_bc_inflow) % inflow side
    Mobj.sponge_fac(Mobj.nSponge) = SpongeCoeffInflow;
    Mobj.sponge_rad(Mobj.nSponge,1:nObcs) = SpongeRadInflow;
  else
    Mobj.sponge_fac(Mobj.nSponge) = SpongeCoeffOutflow;
    Mobj.sponge_rad(Mobj.nSponge,1:nObcs) = SpongeRadOutflow;
  end;
  
end;
write_FVCOM_sponge(Mobj,spgfile);


% write mhke file
write_FVCOM_mhke(mhke_x,mhke_y,mhke_At,mhke_Cp,mhkfile);

% report estimate minimum time step
Mobj = estimate_ts(Mobj,2.5,1);
fprintf('estimated minimum time step in seconds: %f\n',min(Mobj.ts));

% plot everything
%plot_field(Mobj,Mobj.h,'title','domain','withextra',true,'showgrid','false');

%dump_dtascii([casename modelid],Mobj);

% dump node lists for plotting vectors on regular grids
if(make_node_lists);
  
  for ll=1:numel(vecres);
  dx = vecres(ll);
  [X,Y] = meshgrid(min(Mobj.x):dx:max(Mobj.x),min(Mobj.y):dx:max(Mobj.y));
  [ni,nj] = size(X);
  cnt = 0;
  pt = zeros(ni*nj,1);
  for j=1:nj
  for i=1:ni
     cnt = cnt + 1;
     xpt = X(i,j);
     ypt = Y(i,j);
     [mind,pt(cnt)] = min( sqrt( (Mobj.x-xpt).^2 + (Mobj.y-ypt).^2));
  end;
  end;
  pt = unique(pt);
  cnt = numel(pt);
  fid = fopen(['wh' modelid '_vecnodes' num2str(dx) '.dtascii'],'w');
  fprintf(fid,'%s\n','Seq_Vector Nodes');
  fprintf(fid,'%s\n','String');
  fprintf(fid,'%s\n','NumberList');
  fprintf(fid,'%s\n','Vector Nodes');
  fprintf(fid,'%s\n','Real');
  fprintf(fid,'%d %d %d\n',cnt,1,1);
  for i=1:cnt
    fprintf(fid,'%d\n',pt(i));
  end;
  fclose(fid);
  end;
end;
fprintf('mesh has %d vertices \n',Mobj.nVerts);


