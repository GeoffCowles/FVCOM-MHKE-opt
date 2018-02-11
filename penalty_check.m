function Pen = penalty_check(del_x,del_y)

%del_x= 231.8619;
%del_y=21.2633;
optparams

if del_y >= lat_spacing_total
    Pen=0;
    
else
del_x_actual= ((lat_spacing_total-del_y)*long_spacing_total)/lat_spacing_total;

if del_x < del_x_actual
   Pen= penaltymin(del_x,del_x_actual,del_x_actual);
else
    Pen= 0;
end
end
