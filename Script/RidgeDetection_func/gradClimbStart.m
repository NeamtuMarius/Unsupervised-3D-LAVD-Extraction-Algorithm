function [x0,y0,z0] = gradClimbStart(thr0EigPrc,thrdVmPrc,thrVmPrc,Vm,xm,ym,zm)
% gradClimbStart helps in defining the startig position for gradient climbing 
% threshold evaluation is defined as a percentile of the values of the LAVD
% field

% writing features
xv = unique(xm); yv = unique(ym); zv = unique(zm);
dx = mean(diff(xv)); dy = mean(diff(yv)); dz = mean(diff(zv));

% evaluating gradient & hessian & its eigs
[Vmx Vmy Vmz] = gradient(Vm,dx,dy,dz); % gradient
[Vmxx Vmxy Vmxz] = gradient(Vmx,dx,dy,dz); % laplacian
[Vmyx Vmyy Vmyz] = gradient(Vmy,dx,dy,dz);
[Vmzx Vmzy Vmzz] = gradient(Vmz,dx,dy,dz);
dVm = sqrt(Vmx.^2+Vmy.^2+Vmz.^2);
[e1,e2,e3] = EigEval3D(Vmxx,Vmxy,Vmxz,Vmyx,Vmyy,Vmyz,Vmzx,Vmzy,Vmzz);

% finding thresholds
thr0Eig = prctile(abs(e1(:)),thr0EigPrc);
thrdVm = prctile(dVm(:),thrdVmPrc);
thrVm = prctile(Vm(:),thrVmPrc);

% defining initial position:
% finding point with: i) 1st eigval ~ 0; ii) 2nd and 3rd eigvals < 0; iii) high value of Vm
% ind = e2(:) < -thr0Eig  & e3(:) < -thr0Eig & Vm(:)>thrVm;
ind = abs(e1(:)) < thr0Eig & e2(:) < 0 & e3(:) < 0 & Vm(:)>thrVm & dVm(:)>thrdVm;
% ind = e2(:) < -thr0Eig  & e3(:) < -thr0Eig & dVm(:)>thrdVm;
% ind = e2(:) < -thr0Eig  & e3(:) < -thr0Eig & abs(e1(:)) < thr0Eig;
x0 = xm(ind);
y0 = ym(ind);
z0 = zm(ind);

end

