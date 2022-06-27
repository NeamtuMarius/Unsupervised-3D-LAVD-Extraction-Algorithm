function [x0,y0,z0] = gradClimbStartOld(thr0EigPrc,thrdVmPrc,thrVmPrc,e1,e2,e3,dVm,Vm,xm,ym,zm,cpFlag)
%gradClimbStart help in defining the startig position for gradient climbing 
%   Detailed explanation goes here
% threshold evaluation as percentile
thr0Eig = prctile(abs(e1(:)),thr0EigPrc);
thrdVm = prctile(dVm(:),thrdVmPrc);
thrVm = prctile(Vm(:),thrVmPrc);

% evaluating the gradient




% defining initial position:
% finding point with: i) 1st eigval ~ 0; ii) 2nd and 3rd eigvals < 0; iii) high value of Vm
ind = e2(:) < -thr0Eig  & e3(:) < -thr0Eig & Vm(:)>thrVm;
% ind = abs(e1(:)) < thr0Eig & e2(:) < -thr0Eig  & e3(:) < -thr0Eig & Vm(:)>thrVm;
% ind = e2(:) < -thr0Eig  & e3(:) < -thr0Eig & dVm(:)>thrdVm;
% ind = e2(:) < -thr0Eig  & e3(:) < -thr0Eig & abs(e1(:)) < thr0Eig;
x0 = xm(ind);
y0 = ym(ind);
z0 = zm(ind);


if cpFlag == 1
figure % control figure for gradient ascend
p = patch(isosurface(xm,ym,zm,Vm,thrVm));
isonormals(xm,ym,zm,Vm,p); hold on
set(p,'facecolor','red','edgecolor','none','FaceAlpha',0.5);
camlight; lighting gouraud;
plot3(x0(:),y0(:),(z0(:)),'k.');
xlim([min(xm(:)) max(xm(:))]);
zlim([min(zm(:)) max(zm(:))]);
ylim([min(ym(:)) max(ym(:))]);
camproj('perspective')
title('control plot for gradient ascend starting positions')
legend('isoVm','starting pos','Location','SE')
grid off
box on
axis equal
view(3)
end


end

