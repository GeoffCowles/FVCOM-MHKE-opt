function [tpow,powF,powE] = run_case_mesh(mhke_x,mhke_y)


% set optimization parameters
optparams

% remove the DONE file if it exists
if(hydra_run)
  if(exist('DONE'))
    delete('DONE')
  end;
end;


if(hydra_run)
  run_command = 'qsub runscript1';
  ncfile = which('ic01_0001.nc');
 else
  run_command = '/opt/local/bin/mpiexec -np 1 ./fvcom --casename=ic01';
  ncfile = which('ic01_0001.nc');
end

%-------------------------------------------------------------------
% Flood
%-------------------------------------------------------------------

% generate the setup
make_idchan_model(mhke_x,mhke_y,1);


% execute the run
t1 = tic;
disp('running FVCOM flood simulation');
pause(5);
system(run_command);
disp(run_command);

if(hydra_run)
  cycle_until_done;
end;
fprintf('FVCOM execution required %d seconds \n',ceil(toc(t1)));


% read the power from the output
f = dir(ncfile);
if(f.bytes < 1e6)
  disp('something is wrong with FVCOM sim');
  disp('no power computed');
  powF = [-99];
else
  time = ncread(ncfile,'time');
  if(time(end) ~= Tend)
    disp('FVCOM run exited early')
    powF = -99*ones(numel(time),1); 
  else
    powF = ncread(ncfile,'MHKE_POWER'); 
  end;
  ntimes = numel(time);
  %if(hydra_run)
    %system(['/usr/local/bin/ncks -FO -d time,' num2str(ntimes) ',' num2str(ntimes) ',1 ic01_0001.nc ic01_' num2str(caseid) 'flood.nc' ]);
  %end;
end;
tpow = powF(end);

if(one_way); return; end;

%-------------------------------------------------------------------
% Ebb 
%-------------------------------------------------------------------

% generate the setup
make_idchan_model(mhke_x,mhke_y,0);


% execute the run
t1 = tic;
disp('running FVCOM ebb simulation');
pause(5);
system(run_command); 
if(hydra_run)
  cycle_until_done;
end;
fprintf('FVCOM execution required %d seconds \n',ceil(toc(t1)));

% read the power from the output 
f = dir(ncfile); 
if(f.bytes < 1e6)
  disp('something is wrong with FVCOM sim');
  disp('no power computed');
  powE = [-99];
else
  time = ncread(ncfile,'time'); ntimes = numel(time);
  if(time(end) ~= Tend)
    disp('FVCOM run exited early')
    powE = -99*ones(numel(time),1); 
  else
    powE = ncread(ncfile,'MHKE_POWER'); 
  end;
  %if(hydra_run);
    %system(['/usr/local/bin/ncks -FO -d time,' num2str(ntimes) ',' num2str(ntimes) ',1 ic01_0001.nc ic01_' num2str(caseid) 'ebb.nc' ]);
  %end;
end;
tpow = 0.5*(powF(end)+powE(end));
%plot(time,powF,'r'); hold on; plot(time,powE,'k');

