%% main_EulerianLAVDCalculator
% This main computes the ridges of the LAVD field based on a Gradient
% Climbing Algorithm

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


%% gradient ascend starting point
% percentage treshold for defining initial position in gradient ascend.
thr0EigPrc = 10;    % must be close to 0 <--
thrdVmPrc = 80;     % threshold on Vm gradient. 
thrVmPrc = 50;      % threshold on Vm.

[x0,y0,z0] = gradClimbStart(thr0EigPrc,thrdVmPrc,thrVmPrc,Vm,xm,ym,zm);


%% gradient ascend
stepNum = 6000;     % number of maximum time steps
relDiff = 0.075;      % Relative difference between two subsequent binary images in percentage. 

[xpp ypp zpp] = gradAsc(x0,y0,z0,stepNum,xm,ym,zm,Vm,relDiff,1);
points = [xpp ypp zpp];
points(any(isnan(points),2),:) = []; % remove nan values


% control figure
if figFlag ==1
figure('Renderer','opengl','units','normalized','outerposition',[0 0 1 1])
plot3(points(:,1),points(:,2),points(:,3),'k.','MarkerSize',2); hold on
daspect([1 1 1])
box on
end


%% clusterize points 
distmin = mean(diff(unique(xm)))*4; % minimum distance between points to belong to the same cluster
thrLength = 30; % minimum number of points to define a cluster

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
dsmin = 1*(sqrt(mean(diff(unique(xm))).^2 + mean(diff(unique(ym))).^2 + mean(diff(unique(zm))).^2));
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
minlr = 5; % *minimum number of left and right points to discard
ws = 30; % *window size for the spline interpolation

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
ridges.lines  = linesf;
ridges.thr0EigPrc = thr0EigPrc;
ridges.thrVmPrc = thrVmPrc;
ridges.thrdVmPrc = thrdVmPrc;
ridges.thrLength = thrLength;
ridges.dsmin = dsmin;

if saveFlag ==1; save(['./data/Ridges/RidgesGradientClimbing_' expName '.mat'],'ridges'); end









