function new_par = boundary_constraint(pos)

% set the parameters controlling the run
optparams

xvals_boundary = [pos(1:2:end)];
yvals_boundary= [pos(2:2:end)];

global yturb;

varhi_new_y= yturb;
varlo_new_y= -yturb;

global npar;
new_par= zeros(1,npar);
%%
%Newly added boundary contraint
dummy= npar/2;

for c= 1:dummy
   
    if (yvals_boundary(c)>=yturb || yvals_boundary (c)<=-yturb)
        
      yvals_boundary(c)= (varhi_new_y-varlo_new_y)*rand+varlo_new_y;
       
     
    end
    p=(2*c)-1;
    new_par(1,p)= xvals_boundary(c);
    q= 2*c;
    new_par(1,q)= yvals_boundary(c);
    end
