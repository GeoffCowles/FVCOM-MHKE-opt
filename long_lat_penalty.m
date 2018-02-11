function [distP,penalty_overall,penalty_violate] = long_lat_penalty(x,y)
%%
optparams

%xhi=878500;   %design variable upper limit
%xlo=878000;   %design variable lower limit
%yhi=-167850;
%ylo=-166700;

%laxc=[-8e3 8e3];
%layc=[0 0];
%x=(xhi-xlo)*rand(1,5)+xlo;
%y=(yhi-ylo)*rand(1,5)+ylo;
%%
%x=[-10.65995 -50.06609 -96.52533 53.16803 132.4482 -37.33607 -366.6381 -907.9775 149.8764 -600.2939];
%y=[20.80022 -41.08257 -248.0341 -179.0862 73.89728 114.1782 -110.4295 -321.7709 -132.5573 311.1172];


                                                                                                  
%scatter(x,y)
%hold on

%line (laxc,layc)
%hold on

v = [laxc(2)-laxc(1) layc(2)-layc(1)];
t= [0 1];
a= numel(t);
lpt_x=zeros(a,1);
lpt_y=zeros(a,1);
rpt_x=zeros(a,1);
rpt_y=zeros(a,1);
%%
for i=1:a
   
  lpt_x(i)= laxc(i)-v(2);
  lpt_y(i)= layc(i)+v(1);
  rpt_x(i)= laxc(i)+v(2);
  rpt_y(i)= layc(i)-v(1);
  %xlines= [lpt_x(i) rpt_x(i)];
  %ylines= [lpt_y(i) rpt_y(i)];
  %line(xlines,ylines)
end
%%
slope_long_axis= (layc(2)-layc(1))/(laxc(2)-laxc(1));
slope_lat_axis= (rpt_y(2)-lpt_y(2))/(rpt_x(2)-lpt_x(2));
aa= atan(slope_long_axis);
bb= atan(slope_lat_axis);


long_axis_new_x=[];
long_axis_new_y=[];
lat_axis_new_x=[];
lat_axis_new_y=[];
long_axis_new_x_mirror=[];
long_axis_new_y_mirror=[];
lat_axis_new_x_mirror=[];
lat_axis_new_y_mirror=[];
xx=[];
yy=[];
z=numel(x);
global penalty_overall
global penalty_violate

penalty_overall=zeros(1,z);
%%
%j=2;
%%
for j=1:z
    %%
    long_axis_new_x(j)= (long_spacing_total*cos(aa))+x(j);
    long_axis_new_y(j)= (long_spacing_total*sin(aa))+y(j);
    lat_axis_new_x(j)= (lat_spacing_total*cos(bb))+x(j);
    lat_axis_new_y(j)= (lat_spacing_total*sin(bb))+y(j);
    
    %xl= [x(j) long_axis_new_x(j)];
    %yl= [y(j) long_axis_new_y(j)];
    %line(xl,yl)
    
    %xl_2= [x(j) lat_axis_new_x(j)];
    %yl_2= [y(j) lat_axis_new_y(j)];
    %line(xl_2,yl_2)
    %%
   long_axis_new_x_mirror(j)= x(j)-(long_spacing_total*cos(aa));
   long_axis_new_y_mirror(j)= y(j)-(long_spacing_total*sin(aa));
   lat_axis_new_x_mirror(j)= x(j)-(lat_spacing_total*cos(bb));
   lat_axis_new_y_mirror(j)= y(j)-(lat_spacing_total*sin(bb));

%xlm= [x(j) long_axis_new_x_mirror(j)];
%ylm= [y(j) long_axis_new_y_mirror(j)];
%line(xlm,ylm)
%xlm_2= [x(j) lat_axis_new_x_mirror(j)];
%ylm_2= [y(j) lat_axis_new_y_mirror(j)];
%line(xlm_2,ylm_2)

%%
   long_axis_overall_x=[long_axis_new_x(j) long_axis_new_x_mirror(j)];
   long_axis_overall_y=[long_axis_new_y(j) long_axis_new_y_mirror(j)];

%long_axis_overall_x=[long_axis_new_x_mirror(j) long_axis_new_x(j)];
%long_axis_overall_y=[long_axis_new_y_mirror(j) long_axis_new_y(j)];

lat_axis_overall_x= [lat_axis_new_x(j) lat_axis_new_x_mirror(j)];
lat_axis_overall_y= [lat_axis_new_y(j) lat_axis_new_y_mirror(j)];

%lat_axis_overall_x= [lat_axis_new_x_mirror(j) lat_axis_new_x(j)];
%lat_axis_overall_y= [lat_axis_new_y_mirror(j) lat_axis_new_y(j)];
%%

xv=[long_axis_new_x(j) lat_axis_new_x_mirror(j)...
        long_axis_new_x_mirror(j) lat_axis_new_x(j)];
yv= [long_axis_new_y(j) lat_axis_new_y_mirror(j)...
        long_axis_new_y_mirror(j) lat_axis_new_y(j)];
    
%fill(xv,yv,'w')
%%
u= x;
u(:,j)=[];
xq=u;
v= y;
v(:,j)=[];
yq=v;
[in,on] = inpolygon(xq,yq,xv,yv);
%%
if numel(xq(in))==0
     penalty_overall(j)= 0;     
else
 %%   
 xx=xq(in);
 yy=yq(in);
 scatter(xx,yy,'r');
 %%
 n=numel(xx);
 penalty_violate=[];
   %% 
 for k=1:n
    %%
  l=xx(k);
  m=yy(k);
        
  violate_point= [l m];
  e= long_axis_overall_x(2)-long_axis_overall_x(1);
  g= long_axis_overall_y(1)- violate_point(2);
  h= long_axis_overall_x(1)- violate_point(1);
  q= long_axis_overall_y(2)-long_axis_overall_y(1);
  del_y= abs((e*g)-(h*q))/sqrt(e^2+q^2);
    
  e2=   lat_axis_overall_x(2)-lat_axis_overall_x(1);
  g2= lat_axis_overall_y(1)- violate_point(2);
  h2= lat_axis_overall_x(1)- violate_point(1);
  k2= lat_axis_overall_y(2)-lat_axis_overall_y(1);
  del_x= abs((e2*g2)-(h2*k2))/sqrt(e2^2+k2^2);
    
  penalty_violate(k)=penalty_check(del_x,del_y);    
      
 end
%%
  penalty_overall(j)= sum(penalty_violate);
end
end
distP= sum(penalty_overall);
     
