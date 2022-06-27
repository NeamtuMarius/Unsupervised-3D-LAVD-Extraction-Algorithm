function [ring] = getRing(u,x0,dd,xm,ym,zm,Vm,Vmthr,nddsearch,DeficiencyThresh)
%get3Anulus Summary of this function goes here
%   u and x0 are the unit vector and point defining the plane
%   xm,ym,zm,Vm is the coordinates and field
%   Vmthr is the iso-value
%   OUTPUT:
% 	anulus.center   = aaverage verage of anulus coordinates;
%   anulus.ver      = anulus coordinates;
%   anulus.u        = unit vector defining the plane;
%   anulus.x0       = point defining the plane;
%   anulus.isoVm    = threshold value;
%   anulus.isClose  = flag of closeness;
%   anulus.verPlane = planar anulus coordinates;
%   anulus.isCovex  = flag of convexity;    

% orthonormal basis for bringing the 3d anulus in planan coordinates
u1 = [-u(2), u(1), 0];                          u1 = u1/vecnorm(u1);
u2 = [-u(3)*u(1), -u(2)*u(3), u(1)^2 + u(2)^2]; u2 = u2/vecnorm(u2);
u3 = u;                                         u3 = u3/vecnorm(u3);

% extracting anulus
[x_plane,y_plane,z_plane] = getPlane(u(1),u(2),u(3),x0(1),x0(2),x0(3),dd(1),dd(2),dd(3),nddsearch*min(dd)); % get the plane orth to u and passing per x0
if isnan(x_plane) == 1; blankOutput; return; end % return if it didnt work

% ring = Vm_plane;
s = contourslice(xm,ym,zm,Vm,x_plane,y_plane,z_plane,[Vmthr Vmthr]); % get 3D isolines

if length(s) == 0; blankOutput; return; end % return if it didnt work
% disp([num2str(length(s)) ' ring(s) founded'])

for i = 1:length(s)
    % disp(['processing ' num2str(i) ' ring of ' num2str(length(s))])
    ver = s(i).Vertices; ver(end,:) = [];
    verPlane(:,1) = (ver-x0)*u2';
    verPlane(:,2) = (ver-x0)*u1';
    
    verv{i} = ver;
    isClose(i)  = (sum(ver(1,:) - ver(end,:)) == 0);
    isInside(i) = inpolygon(0,0,verPlane(:,1),verPlane(:,2));
    isConvex(i) = IsContourConvex(verPlane(:,1),verPlane(:,2),DeficiencyThresh);
    
    clear verPlane
    
end

% indices for close and inside points
indGood = find(isInside.*isClose == 1);

if length(indGood) == 0; blankOutput; return; end % return if it didnt work

% structure building
ring.ver = verv{indGood(1)};
ring.isConvex = isConvex(indGood(1));
ring.t = u;
ring.x0 = x0;
close


function blankOutput
    ring.ver = [nan nan nan];
    ring.isConvex = [nan];
    ring.t = u;
    ring.x0 = x0;
    close
    
end


end

