function [omegax omegay omegaz] = gradientCalc(u,v,w,dx,dy,dz,dt) 

% gradients evaluation
[ux uy uz ut] = gradient(u,dx,dy,dz,dt); % velocity gradient tensor
[vx vy vz vt] = gradient(v,dx,dy,dz,dt);
[wx wy wz wt] = gradient(w,dx,dy,dz,dt);

% Sxx = (ux + ux)/2; Sxy = (uy + vx)/2; Sxz = (uz + wx)/2; % simmetric (strain rate tensor)
% Syx = (vx + uy)/2; Syy = (vy + vy)/2; Syz = (vz + wy)/2;
% Szx = (wx + uz)/2; Szy = (wy + vz)/2; Szz = (wz + wz)/2;

omegax = wy - vz; % vorticity
omegay = uz - wx;
omegaz = vx - uy;

% omega2 = omegax.^2 + omegay.^2 + omegaz.^2; % enstrophy


end
