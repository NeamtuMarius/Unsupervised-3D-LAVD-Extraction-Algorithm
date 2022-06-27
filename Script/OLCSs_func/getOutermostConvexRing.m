function [ring] = getOutermostConvexRing(u,x0,dd,xm,ym,zm,Vm,dVmIso,nddsearch,DeficiencyThresh)
%getOutermostConvexRing Summary of this function goes here
%   Detailed explanation goes here
x = unique(xm); y = unique(ym); z = unique(zm);
Vm_interp = griddedInterpolant({x,y,z},permute(Vm,[2,1,3]),'spline','none');
VmIso0 = Vm_interp(x0(1),x0(2),x0(3));

VmIso = VmIso0;
isConvex = 1;
while isConvex == 1
    VmIso = VmIso - dVmIso;
    ring = getRing(u,x0,dd,xm,ym,zm,Vm,VmIso,nddsearch,DeficiencyThresh);
    isConvex = ring.isConvex;
    
end






end

