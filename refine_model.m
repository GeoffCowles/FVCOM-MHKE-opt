
function [Mobj] = refine_model(Mobj,mark,nsmooth,tol)
%==============================================================================
% Refine an FVCOM mesh object
%
%   Input:
%     Mobj:  FVCOM mesh object
%     mark:  Vector of size Mobj.nElems, =1 refine
%     nsmoo: number of smoothing iterations (optional, default=20)
%     tol:   smoothing tolerance (optional, default=.05) 
%  
%   TODO
%     Bathymetry will be interpolated on to the new mesh using the bathy
%     on the old mesh.  However you may wan to reinterpolate from the bathy
%     database in order to get new features.
%
%     Coriolis will be interpolated to new mesh using old
%
%     Open boundary node node list will be adjusted if new points are 
%     inserted in the open boundary. 
%  
%     Caution:  this may affect the numbering of river nodes
%     Caution:  if new point is added on open boundary, open boundary
%               forcing files may need to be adjusted. 
%
%==============================================================================


% refine the mesh
[p,t] = refine([Mobj.x Mobj.y],Mobj.tri,mark);

% smooth the mesh
[p,t] = smoothmesh(p,t,nsmooth,tol);

% fill mesh object
Mobj.x = p(:,1);
Mobj.y = p(:,2);
Mobj.tri = t;

% set dimensions
Mobj.nElems = size(t,1);
Mobj.nVerts = size(p,1);


