addpath(genpath('../dependencies'));
path_to_gmsh = '/usr/local/bin/gmsh';

plotit = 1; %=1 to plot turbine positions, =0 for no graphics
hydra_run = 1; %=1 running on cluster, = 0 running locally

% optimization parameters
nTurb = 5;
one_way = 0; %=1 only consider flood, = 0, average power from flood/ebb
MaxIt=60; %max iterations
popsize=10; %set population size
mutrate= 0.2; %set mutation rate
selection=0.5; %fraction of population kept
plotint = 100;
Nodes=4;

Psep_mag= -1000;
mind=500; 

% set geometry parameters
chan_length   = 5000;
chan_width    = 1000;
constriction  = 0.8; %=1, no constriction, =0.5, channel constricted 50%
fillet_radius = 500;
ob_radius     = 5000;
ds_channel    = 75; %mesh spacing in the channel
ds_outer      = 2000.; %mesh spacing at the open boundary

% set bathymetry parameters 
chan_depth      = 20.; %depth in the channel
max_depth       = 75.; %depth at open boundary
bath_smooth_its = 5;

% sponge damping coefficient and radius
SpongeCoeffInflow = 1e-4;
SpongeCoeffOutflow = 1e-3;
SpongeRadInflow   = 2*ds_outer;
SpongeRadOutflow   = 2*ds_outer;


%MHKE parameters 
At  = 600; %turbine cross sectional area
Cp  = 0.5; %turbine power coefficient
refine_mhke=0; %=1 to create refinement area, = 0 for no refinement
dxR = 1000; %area to refine in x direction on either side of turbine
dyR = 100; %area to refine in y direction on either side of turbine

% head (free surface elevation difference between inflow and outflow in m
head = 0.23;


% parameters related to FVCOM run 
Tend = 0.25;  %end time of the FVCOM simulation in days

