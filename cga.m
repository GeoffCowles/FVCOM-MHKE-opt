%-------------------------------------------------------------------------
%This is the main operating code for continuous genetic algorithm.
% Selection Method Used: Roulette Wheel Selection
% Crossover Method Used: Single Point Crossover
% For Fundamentals of Continuous Genetic Algorithm: Refer to "Book by Haupt"
% in the resource folder (Page 51) 
% Matlab Source Code: Refer to "carrGenet" in the resource folder (Page 33) 
%-------------------------------------------------------------------------
%%
close all;
clear all;
clc;
clear;
tic
addpath('../Mesh2Dv24');

% set the parameters controlling the run
optparams

check_init = 0;
%figure('units','normalized','outerposition',[0 0 0.5 0.5])


global npar;
npar=nTurb*2; %number of optimization variables

dy = constriction*chan_width/npar;

global yturb; 
yturb = 0.5*constriction*chan_width-dy/2;



varhi=chan_length*.5;   %design variable upper limit
varlo=-chan_length*.5;   %design variable lower limit

%Stopping criteria
maxit=MaxIt;  %max number of iterations


Nt = npar;

%population members that survive in each iterations
keep=floor(selection*popsize);

%total number of mutations that occurs in each iterations
nmut=ceil ((popsize-1)*npar*mutrate);

M=ceil((popsize-keep)/2); %number of pairs to mate in each iterations

% Create the initial population
iga=0; %generation counter initialized

% Randomly generating initial population
par= (varhi-varlo)*rand(popsize,npar)+varlo;

for n=1:popsize
  pos= par(n,:);
  new_par = boundary_constraint(pos);
  par(n,:)= new_par;
end

 
count=0;

[cost]=zeros(1, popsize);


for s= 1:Nodes
      ff= 'Node';
      folder = strcat(ff, num2str(s)); 
      mkdir('..\',folder);
      copyfile ('..\Ideal_Channel_Test_Cluster',['..\',folder]);
end

if matlabpool('size') > 0 % checking to see if my pool is already open
    matlabpool close
end

%if ~strcmp(getenv('PBS_JOBID'),'')
 % sched = findResource('scheduler','type','local');
 % local_scheduler_data=[sched.DataLocation,'/',getenv('PBS_JOBID')]
 % mkdir(local_scheduler_data);
 % sched.DataLocation=local_scheduler_data;
%end

matlabpool (Nodes);


spmd
    cd ../
    cd(sprintf('Node%d', labindex));
end

% Get cost for initial population
parfor (f=1:popsize,Nodes)
  pos= par(f,:);
  xvals= pos(1:2:end);
  yvals= pos(2:2:end);
  cost(f) = 1*run_case(xvals,yvals);
  count=count+1;
end;

matlabpool close


% Max cost in element 1 & sorting cost values in descending order
[cost,ind]=sort(cost,'descend'); 

par=par(ind,:); %sorting population based on sorted cost values

% For Post Processing Plotting
best= par(1,:);
best_xvals= best(1:2:end);
best_yvals= best(2:2:end);
best_cost= 1*run_case_mesh(best_xvals,best_yvals);
count=count+1;

global Mobj;
mesh{1}= Mobj;
Coords{1}=par; % Keeping trace of population in each iteration
Cost_val{1}= cost;
k(1)= count;

maxc(1)=max(cost); % maxc contains max of population for reporting
meanc(1)=mean(cost);  % meanc contains mean of population for reporting
time(1)=toc;
%%

%-------------------------------------------------------------------------
%Iterate through generations (Main Loop)
%------------------------------------------------------------------------
while iga<maxit
  iga=iga+1;  %increment the generation counter
  
  %---------------------------------------------------------------
  %Selection of Mother & Father by Roulette Wheel Selection Method
  %For Fundamentals: Refer to "Roulette-wheel selection" in resource folder
  %---------------------------------------------------------------
  
  prob=flipud([1:keep]'/sum([1:keep])); % weights chromosomes
  odds=[0 cumsum(prob(1:keep))']; % probability    distribution   function
  
  pick1=rand(1,M); % mate 1(vector of length M with randoms between 0 & 1)
  pick2=rand(1,M); % mate 2
  
  %ma and pa contain the indices of the chromosomes that will mate and ...
  % choosing  integer k with  probability p(k)
  
  ic=1;
 
while ic <=M
    
  for id=2:keep+1
       
     if pick1(ic)<=odds(id) && pick1(ic)>odds(id-1)
        ma(ic)=id-1;
     end
     
     if pick2(ic)<=odds(id) && pick2(ic)>odds(id-1)
        pa(ic)=id-1;
     end
  end
ic=ic+1;
end

%-------------------------------------------------------------------------
%End of Selection Method
%------------------------------------------------------------------------


%-------------------------------------------------------------------------
%Performs mating using single point crossover
%------------------------------------------------------------------------
% Indexing needed to address offsprings

ix=[];
for cc=1:M
    
ix(cc)=(2*cc)-1;
end
  

xp=ceil(rand(1,M)*(npar)); % Randomly generated crossover point
r=rand(1,M);  % Mixing parameter used with crossover point element

for ic=1:M
  
  % Mating of ma and pa at crossover point
  xy=par(ma(ic),xp(ic))-par(pa(ic),xp(ic)); 
    
  par(keep+ix(ic),:)=par(ma(ic),:); %1st offspring from ma
  par(keep+ix(ic)+1,:)=par(pa(ic),:); %2nd offspring from pa
  
  %Blending of mixing parameter and mating parameter in 1st offspring
  par(keep+ix(ic),xp(ic))=par(ma(ic),xp(ic))-r(ic).*xy;
  
  %Blending of mixing parameter and mating parameter in 2nd offspring
  par(keep+ix(ic)+1,xp(ic))=par(pa(ic),xp(ic))+r(ic).*xy;
  
  %Swapping the remaining elements after the crossover point
  if xp(ic)<(npar)
      
    %1st offspring
    par(keep+ix(ic),:)=[par(keep+ix(ic),1:xp(ic)) par(keep+ix(ic)+1,xp(ic)+1:(npar))];
    
    %2nd offspring
    par(keep+ix(ic)+1,:)=[par(keep+ix(ic)+1,1:xp(ic)) par(ma(ic),xp(ic)+1:(npar))];
    %par(keep+ix(ic)+1,:)=[par(keep+ix(ic)+1,1:xp(ic)) par(keep+ix(ic),xp(ic)+1:(npar))];
  
  end
end

par=par(1:popsize,:);

%-------------------------------------------------------------------------
%End of Crossover Method
%------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Mutate the population
%------------------------------------------------------------------------

% Selection of rows for mutation
mrow=sort(ceil(rand(1,nmut)*(popsize-1))+1);

% Selection of rows for mutation
mcol=ceil(rand(1,nmut)*(Nt));
  
for ii=1:nmut
 
 % Picking randomly generated number at a particular row and column
 % selected for mutation
 par(mrow(ii),mcol(ii))=(varhi-varlo)*rand+varlo;  
 
end

%-------------------------------------------------------------------------
%End of Mutation
%------------------------------------------------------------------------

%The new offspring and mutated chromosomes are evaluated

for n=1:popsize
  pos= par(n,:);
  new_par = boundary_constraint(pos);
  par(n,:)= new_par;
end

[cost]=zeros(1, popsize);


if matlabpool('size') > 0 % checking to see if my pool is already open
    matlabpool close
end

if ~strcmp(getenv('PBS_JOBID'),'')
  sched = findResource('scheduler','type','local');
  local_scheduler_data=[sched.DataLocation,'/',getenv('PBS_JOBID')]
  mkdir(local_scheduler_data);
  sched.DataLocation=local_scheduler_data;
end

matlabpool (Nodes);

spmd
    
    cd ../
    cd(sprintf('Node%d', labindex));
   
end    

% Get cost for initial population
parfor (f=1:popsize,Nodes)
  pos= par(f,:);
  xvals= pos(1:2:end);
  yvals= pos(2:2:end);
  cost(f) = 1*run_case(xvals,yvals);
  count=count+1;
end;

matlabpool close

%Sorting the costs and associated population in a descending order
  
[cost,ind]=sort(cost,'descend');
par=par(ind,:);

%Need to get the mesh file from the best offspring
best= par(1,:);
best_xvals= best(1:2:end);
best_yvals= best(2:2:end);
best_cost= 1*run_case_mesh(best_xvals,best_yvals);
count=count+1;

global Mobj;
mesh{iga+1}= Mobj;
Coords{iga+1}=par;
Cost_val{iga+1}= cost;
  
%Do statistics for a single nonaveraging run
maxc(iga+1)=max(cost);
meanc(iga+1)=mean(cost);
k(iga+1)= count;

save('pop.mat','Coords');
save('cost.mat','Cost_val');
save('max_cost.mat','maxc');
save('mesh.mat','mesh');
save('func_eval.mat','k');

  
%Stopping criteria
  
if iga>maxit
    break
end
[iga cost(1)];

%-------------------------------------------------------------------------
% Plot positions in every iterations
%------------------------------------------------------------------------

clf

% For Flood Function
global Mobj;
 patch('Vertices',[Mobj.x,Mobj.y],'Faces',Mobj.tri,...
     'edgecolor','green','FaceColor','none');
      colorbar;
  hold on;
  caxis([1.5,2.5]);
  axis([-8000,8000,-5000,5000]);
  


[maxcost,imax] = max(cost);
  
    
  for i=1:popsize;
    plot(par(i,1:2:end),par(i,2:2:end),'r+')
  end;
  
plot(par(imax,1:2:end),par(imax,2:2:end),'kd',...
      'MarkerFaceColor','w','MarkerSize',12);

  
fprintf('Generation %d best overall is %f\n',iga,max(cost));
%drawnow;
%pause(0.5);
time(iga+1)= toc;
save('elapsed_time.mat','time');
  
end %iga
%--------------------------------------------
% End Main Loop over Generations
%--------------------------------------------


%Displays Number of Generations Vs Maximum Cost
figure (2)
day=clock;
disp(datestr(datenum(day(1),day(2),day(3),day(4),day(5),day(6)),0))
format  short g
disp(['popsize=' num2str(popsize) ' mutrate=' num2str(mutrate) ' # par=' num2str(npar)])
disp(['#generations='   num2str(iga)   '   best   cost='  num2str(cost(1))])
disp('best solution')
disp(num2str(par(1,:)))
disp('continuous genetic algorithm')
iters=0:length(maxc)-1;
plot(iters,maxc,iters,meanc,'-');
xlabel('generation');
ylabel('cost');
hold on

%Displays Number of Functional Evaluation Vs Maximum Cost
%k=linspace(1,count,MaxIt+1);
figure (3)
plot(k,maxc);
xlabel('Function Evaluation');
ylabel('Maximum Cost');
toc

