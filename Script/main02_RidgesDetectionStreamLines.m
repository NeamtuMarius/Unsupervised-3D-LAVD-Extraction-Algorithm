%% housekeeping
clear all
close all
clc


%% addpath
addpath('./RidgeDetection_func');


%% setting input
expName = 'JHISO';
figFlag = 1;
saveFlag = 1;


%% loading data
load(['./data/LAVDs/LAVD_' expName '.mat'])

Vm  = LAVD.VMatrix;
xm = LAVD.xpi;
ym = LAVD.ypi;
zm = LAVD.zpi;

dx = mean(diff(unique(xm)));
dy = mean(diff(unique(ym)));
dz = mean(diff(unique(zm)));

[Vmx Vmy Vmz] = gradient(Vm,dx,dy,dz); % defining gradient of Vm


%% gradient ascend with streamlines starting point
xm0 = xm(1:2:end,1:2:end,1:2:end); % *defining initial streamlines locations. Here a streamline is released every two grid points (1:2:end)
ym0 = ym(1:2:end,1:2:end,1:2:end);
zm0 = zm(1:2:end,1:2:end,1:2:end);
stepSize = 0.1; % *stepsize in streamlines calculation 
stepNum = 200; % *max number of iterations in stream lines calculation 

verts = stream3(xm,ym,zm,Vmx,Vmy,Vmz,xm0,ym0,zm0,[stepSize stepNum]); % calculating the streamlines
verts = cell2mat(verts'); % convert cell into matrix


%% bin counting of the streamlines points
dd = mean(diff(unique(xm(:))))/5; % *bin size to count streamline points (< flow field res)

xEdges = min(xm(:))-dd/2:dd:max(xm(:))+dd/2; % xedges for binning
yEdges = min(ym(:))-dd/2:dd:max(ym(:))+dd/2; % yedges for binning  
zEdges = min(zm(:))-dd/2:dd:max(zm(:))+dd/2; % zedges for binning
xCenters = xEdges+dd/2;
yCenters = yEdges+dd/2;
zCenters = zEdges+dd/2;
N = histcnd(verts(:,1),verts(:,2),verts(:,3),xEdges,yEdges,zEdges); % calculating 3d binning matrix


% % control figure  % cdf of hystogram counting
% if figFlag ==1
% figure('Renderer','opengl')
% [fN,xN] = ecdf(N(:));
% plot(xN,100*fN,'-'); hold on
% set(gca,'xscale','log')
% xlabel('number of points per bin')
% ylabel('ecdf %')
% end


%% get points
Nprc = 99.99; % *threshold in percentile on the number of points in a bin to define the cluster location (should be close to 100)

Nthr = prctile(N(:),Nprc)+0.5; % get the corresponding value of the percentage
[xmCenters ymCenters zmCenters] = meshgrid(xCenters,yCenters,zCenters); % define the location of the grid
points = [xmCenters(N>Nthr) ymCenters(N>Nthr) zmCenters(N>Nthr)]; % define the points of the clusters
points(any(isnan(points),2),:) = []; % remove nan values


% control figure
if figFlag ==1
figure('Renderer','opengl','units','normalized','outerposition',[0 0 1 1])
plot3(points(:,1),points(:,2),points(:,3),'k.','MarkerSize',2); hold on

xslice = [];  
yslice = mean(unique(ym));
zslice = [];
s = slice(permute(ym,[2 1 3]),permute(xm,[2 1 3]),permute(zm,[2 1 3]),permute(Vm,[2 1 3]),xslice,yslice,zslice);
set(s,'EdgeColor','None')
% set(gca,'ColorScale','Log')

s = isosurface(ym,xm,zm,Vm,prctile(Vm(:),98))
isonormals(Vm,patch(s,'FaceColor','r','FaceAlpha',0.2,'EdgeColor','none'))
view(3)
camlight
lighting gouraud
daspect([1 1 1])
box on
end


%% clusterize points 
distmin = mean(diff(unique(xm)))*5; % minimum distance between points to belong to the same cluster. A distance of 5 flow field resolution was selected.
thrLength = 50; % minimum number of points to define a cluster

% clusters = points2clusters(points,distmin,thrLength);
pointsColor = pcsegdist(pointCloud(points),distmin,NumClusterPoints=[thrLength,Inf]);
clusters = [points double(pointsColor)]; clusters(~pointsColor,:) = [];


% control figure
if figFlag ==1
figure('Renderer','opengl','units','normalized','outerposition',[0 0 1 1])
scatter3(clusters(:,1),clusters(:,2),clusters(:,3),2,clusters(:,4)); hold on
colormap(flipud('colorcube'))
daspect([1 1 1])
box on
end


%% get lines from clusters
dsmin = 1*(sqrt(mean(diff(unique(xm))).^2 + mean(diff(unique(ym))).^2 + mean(diff(unique(zm))).^2)); % threshold on minimum distance bewteen adjacent points along a line
lines = clusters2lines(clusters,dsmin);


% control figure
if figFlag ==1
figure('Renderer','opengl','units','normalized','outerposition',[0 0 1 1])
for i = 1:length(lines)
    plot3(lines{i}(:,1),lines{i}(:,2),lines{i}(:,3),'k-'); hold on
end
daspect([1 1 1])
box on
end


%% spline fitting
minlr = 5; % *minimum number of left and right points to discharge
ws = 20; % *window size for the spline interpolation

for i = 1:length(lines) % cicle on the lines
    linesf{i} = beatSpline(lines{i},minlr,ws); % apply beat interpolation to all lines
end


% control figure
if figFlag ==1
figure('Renderer','opengl','units','normalized','outerposition',[0 0 1 1])
dataPlot = cell2mat(linesf'); % convert cell to matrix for plot
stepPlot = 5; % *step for visualization of tangent vectors
quiver3(dataPlot(1:stepPlot:end,1),dataPlot(1:stepPlot:end,2),dataPlot(1:stepPlot:end,3),...
        dataPlot(1:stepPlot:end,4),dataPlot(1:stepPlot:end,5),dataPlot(1:stepPlot:end,6),'r-'); hold on
plot3(dataPlot(:,1),dataPlot(:,2),dataPlot(:,3),'k-');
daspect([1 1 1])
box on
end


%% saving data structure
ridges.lines = linesf;
ridges.thr0EigPrc = nan;
ridges.thrVmPrc = nan;
ridges.thrdVmPrc = nan;
ridges.thrLength = thrLength;
ridges.dsmin = dsmin;

if saveFlag ==1; save(['./data/Ridges/RidgesStreamlines_' expName '.mat'],'ridges'); end









