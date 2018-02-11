% use global optimization to optimize power
clear all; close all;
if(exist('DONE')); delete('DONE'); end;
system('\rm -r *.nc');

% set case parameters
optparams

% set downstream dimension 
ndown = 10;
dx = .5*chan_length/ndown;
xt = linspace(-.5*chan_length+dx,.5*chan_length-dx,ndown);

% set upstream dimension 
nacross = 5;
dy = .5*chan_width/nacross;
yt = linspace(-.5*chan_width+dy,.5*chan_width-dy,nacross);

% main loop over simulations
if(plotit);
  rectangle('Position',[-0.5*chan_length,-0.5*chan_width,chan_length,chan_width]); drawnow; hold on;
end;
nsim  = 400;
nturb = 6;
for i=1:nsim
  fprintf('running turbine array %d of %d\n',i,nsim);
  ind = randperm(ndown*nacross,nturb);
  idown   = floor((ind-1)/nacross)+1;
  iacross = mod(ind-1,nacross)+1;
  xturb   = xt(idown);
  yturb   = yt(iacross);
  if(plotit); plot(xturb,yturb,'r+'); drawnow; end;
  mhke(i).x = xturb;
  mhke(i).y = yturb;
  [mhke(i).tpower,mhke(i).powF,mhke(i).powE] = run_case(xturb,yturb,i);
  fprintf('max power %f MW\n',mhke(i).tpower*1e-6);
  save mhke mhke
end;
save mhke mhke

  

