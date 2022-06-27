function [IsoConv] = rings2IsoConv(rings)
%rings2extSurf connect the rings to a single surface
%   

IsoConv = struct('vertices',[],'faces',[]);
for i = 1:length(rings)-1
    ver1 = rings{i}.ver;
    ver2 = rings{i+1}.ver;
    u1 = rings{i}.t;
    u2 = rings{i+1}.t;
    
    if isnan(sum(ver1(:)))==0 & isnan(sum(ver2(:)))==0 & sum(u1) ~= 0 & sum(u2) ~= 0
        [IsoConvI] = getLatSurf(ver1(:,1),ver1(:,2),ver1(:,3),ver2(:,1),ver2(:,2),ver2(:,3),u1(:,1),u1(:,2),u1(:,3),u2(:,1),u2(:,2),u2(:,3));
        IsoConv.faces = [IsoConv.faces; IsoConvI.faces+size(IsoConv.vertices,1)];
        IsoConv.vertices = [IsoConv.vertices; IsoConvI.vertices];

    end
    
end


end

