function [IsoConv] = getLatSurf(X1,Y1,Z1,X2,Y2,Z2,nx1,ny1,nz1,nx2,ny2,nz2)

curve=convhull(cat(1,X1,X2),cat(1,Y1,Y2),cat(1,Z1,Z2));
IsoConv.vertices=[cat(1,X1,X2) cat(1,Y1,Y2) cat(1,Z1,Z2)];IsoConv.faces=curve;


aa = IsoConv.vertices(IsoConv.faces(:, 2), :) - IsoConv.vertices(IsoConv.faces(:, 1), :);
bb = IsoConv.vertices(IsoConv.faces(:, 3), :) - IsoConv.vertices(IsoConv.faces(:, 1), :);
cc = cross(aa, bb, 2);
% %%%%% Normals
normals=[];
for jj=1:3 %number of dimensions 
  normals(:,jj)=abs(cc(:,jj)./sqrt(sum(cc.^2, 2)));
end


% PCA1=pca([X1 Y1 Z1]);
% e11=PCA1(:,1)'; e21=PCA1(:,2)' ;e31=PCA1(:,3)';  % 3 principal vector(3 eigenvector) of "input data"
% PCA2=pca([X2 Y2 Z2]);
% e12=PCA2(:,1)'; e22=PCA2(:,2)' ;e32=PCA2(:,3)';  % 3 principal vector(3 eigenvector) of "input data"


tol=1e-2;

% numones=(ismembertol(normals(:,1),abs(e31(1)),tol) & ismembertol(normals(:,2),abs(e31(2)),tol) & ...
%     ismembertol(normals(:,3),abs(e31(3)),tol)) |...
%     (ismembertol(normals(:,1),abs(e21(1)),tol) & ismembertol(normals(:,2),abs(e21(2)),tol) & ...
%     ismembertol(normals(:,3),abs(e21(3)),tol)) |...
%     (ismembertol(normals(:,1),abs(e11(1)),tol) & ismembertol(normals(:,2),abs(e11(2)),tol) & ...
%     ismembertol(normals(:,3),abs(e11(3)),tol)) |...
%     (ismembertol(normals(:,1),abs(e32(1)),tol) & ismembertol(normals(:,2),abs(e32(2)),tol) & ...
%     ismembertol(normals(:,3),abs(e32(3)),tol)) |...
%         (ismembertol(normals(:,1),abs(e22(1)),tol) & ismembertol(normals(:,2),abs(e22(2)),tol) & ...
%     ismembertol(normals(:,3),abs(e22(3)),tol)) |...
%         (ismembertol(normals(:,1),abs(e12(1)),tol) & ismembertol(normals(:,2),abs(e12(2)),tol) & ...
%     ismembertol(normals(:,3),abs(e12(3)),tol));

numones=(ismembertol(normals(:,1),abs(nx1),tol) & ismembertol(normals(:,2),abs(ny1),tol) & ...
    ismembertol(normals(:,3),abs(nz1),tol)) |...
    (ismembertol(normals(:,1),abs(nx2),tol) & ismembertol(normals(:,2),abs(ny2),tol) & ...
    ismembertol(normals(:,3),abs(nz2),tol));

IsoConv.faces=IsoConv.faces(~numones,:);



end
% 
% X1=c_str.xc{1};Y1=c_str.yc{1};Z1=c_str.zc{1};
% X2=c_str.xc{2};Y2=c_str.yc{2};Z2=c_str.zc{2};
%   figure
%   plot3(X1,Y1,Z1),hold on;
%   plot3(X2,Y2,Z2),hold on;
%   axis equal
%   
% figure
% p = patch(IsoConv);
% view(3);
% set(p,'FaceColor',[0.2,0.8,0.2],'EdgeColor','none','FaceAlpha',0.5);
% camlight; lightangle(-94,12); lighting phong
% set(gca,'color',[151,232,255]/255,'fontsize',12); 
% xlabel('X [m]'); ylabel('Y [m]');  zlabel('Z [m]');
% axis tight; daspect([1,1,1]); view(3); hold on;
% 
% IsoConvsmooth=smoothpatch(IsoConv,1,1,1,40000000000)
% 
% IsoConv_s.vertices=SurfaceSmooth(IsoConv.vertices,IsoConv.faces,'VoxSize',.005);
% 
% 
% figure
% p = patch(IsoConvsmooth);
% view(3);
% set(p,'FaceColor',[0.2,0.8,0.2],'EdgeColor','none','FaceAlpha',0.5);
% camlight; lightangle(-94,12); lighting phong
% set(gca,'color',[151,232,255]/255,'fontsize',12); 
% xlabel('X [m]'); ylabel('Y [m]');  zlabel('Z [m]');
% axis tight; daspect([1,1,1]); view(3); hold on;
