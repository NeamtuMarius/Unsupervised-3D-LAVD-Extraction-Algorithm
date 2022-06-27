function [x_plane,y_plane,z_plane] = getPlane(nx,ny,nz,cx,cy,cz,dx,dy,dz,l)

if abs(nx)>abs(ny) && abs(nx)>abs(nz)
    [y_plane,z_plane] = meshgrid((cy-l/2):dy:(cy+l/2),(cz-l/2):dz:(cz+l/2)); 
    x_plane=(-ny.*y_plane-nz.*z_plane+nx.*cx+ny.*cy+nz.*cz)./nx;

elseif abs(ny)>abs(nx) && abs(ny)>abs(nz)
    [x_plane,z_plane] = meshgrid((cx-l/2):dx:(cx+l/2),(cz-l/2):dz:(cz+l/2));
    y_plane=(-nx.*x_plane-nz.*z_plane+nx.*cx+ny.*cy+nz.*cz)./ny;
    
elseif abs(nz)>abs(nx) && abs(nz)>abs(ny)
    [x_plane,y_plane]=meshgrid((cx-l/2):dx:(cx+l/2),(cy-l/2):dy:(cy+l/2));
    z_plane=(-nx.*x_plane-ny.*y_plane+nx.*cx+ny.*cy+nz.*cz)./nz;
    
else
    x_plane = nan;
    y_plane = nan;
    z_plane = nan;
    
end


end

     
      
     
     