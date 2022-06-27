function [OLCS] = getHallerVortex(Vm,xm,ym,zm,x0,t,dd,nddsearch,prcIncr,DeficiencyThresh,ringNumTol)
%getHallerVortex building Haller OLCS 
%   

PERCENTRANK = @(YourArray,TheProbes) reshape(mean(bsxfun(@le,YourArray(:),TheProbes(:).'))*100,size(TheProbes)); % inverse percentage

Vm_interp = griddedInterpolant({unique(xm),unique(ym),unique(zm)},permute(Vm,[2,1,3]),'spline','none');
VmIso0 = min(Vm_interp(x0(:,1),x0(:,2),x0(:,3)));
prcVmIso0 = PERCENTRANK(Vm(:),VmIso0);


VmIso = VmIso0;
prcVmIso = prcVmIso0;

test = 0;
test1 = 0;
isConvex = 1;
while test >= test1 - ringNumTol
    clear isConvex
    prcVmIso = prcVmIso - prcIncr;
    VmIso = prctile(Vm(:),prcVmIso);
    for i = 1:size(x0,1)
        rings{i} = getRing(t(i,:),x0(i,:),dd,xm,ym,zm,Vm,VmIso,nddsearch,DeficiencyThresh);
        isConvex(i,1) = double(rings{i}.isConvex);
        
    end
    
    test = nansum(isConvex);
    test1 = length(isConvex) - sum(isnan(isConvex));
    
    
    
    disp(['at iso ' num2str(VmIso) ' out of ' num2str(VmIso0) ', we have ' num2str(test) ' rings out of '  num2str(test1) ' are convex'])
    
    
end

IsoConv = rings2IsoConv(rings(isConvex == 1));

VmMod = Vm;
VmMod(1,:,:) = 0; VmMod(end,:,:) = 0;
VmMod(:,1,:) = 0; VmMod(:,end,:) = 0;
VmMod(:,:,1) = 0; VmMod(:,:,end) = 0;

ISRF = isosurface(xm,ym,zm,VmMod,VmIso);
ISRF_SP = splitFV(ISRF.faces,ISRF.vertices);

for i = 1:length(ISRF_SP)
    IN = inpolyhedron(ISRF_SP(i),x0,'flipNormals',true);
    control(i) = sum(IN) == size(x0,1);
end
indIN = find(control == 1);
IsoSurf = ISRF_SP(indIN);


% structure building
OLCS.rings = rings;
OLCS.centers = x0;
OLCS.t = t;
OLCS.IsoConv = IsoConv;
OLCS.IsoSurf = IsoSurf;



end

