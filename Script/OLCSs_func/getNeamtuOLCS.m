function [OLCS] = getNeamtuOLCS(Vm,xm,ym,zm,pos,t,dd,nddsearch,percVm0,percVmincr,DeficiencyThresh)
%getNeamtuOLCS building Neamtu OLCS (otermost convex anulus for each point)
%   pos,t of the all ridge points in the curve

for i = 1:size(pos,1)
    x0 = pos(i,:); % point
    u = t(i,:); % unit vector in that point
    
    % anuli{i} = getOutermostConvexAnulus(Vm,xm,ym,zm,x0,u,dd,nddsearch,percVm0,percVmincr,DeficiencyThresh);
    anuli{i} = getOutermostConvexRing(u,x0,dd,xm,ym,zm,Vm,0.5,nddsearch,DeficiencyThresh);
    
    disp([num2str(i/size(pos,1)*100) '%'])

end

% IsoConv = anuli2IsoConv(anuli);

% structure building
OLCS.rings = anuli;
OLCS.centers = pos;
OLCS.directions = t;
% OLCS.IsoConv = IsoConv;

end

