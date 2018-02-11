function [distP] = calc_dist_penalty(x,y,dist)
% compute distance penalty function
% for 1D and 2D case
% G. Cowles, 2016
%


%addpath(pwd);
distP = 0;


    % mesh the turbine positions
    p = [x' y'];
    t = delaunay(x,y);
    [e,~,~,bnd] = connectivity(p,t);
    
    % penalize the edgelengths
    [ne,~] = size(e);
    for i=1:ne
      n1 = e(i,1);
      n2 = e(i,2);
      lngth = sqrt( (x(n1)-x(n2))^2 + (y(n1)-y(n2))^2);
      distP = distP +  penaltymin(lngth,dist,dist);
    end;
    
  
%distP = distP/ne; %should we normalize????
