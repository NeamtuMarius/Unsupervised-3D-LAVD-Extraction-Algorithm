%% main_JHTurbExtractor
% extract data from jht database without any size limitation
% (need to download "turbmat-master" from "https://github.com/idies/turbmat" and add the path of JHT database matlab functions)

%% housekeeping
clear all;
close all;
clc


%% addpath of the JH turb dataset matlab func folder
addpath('./turbmat-master'); % JH functions
addpath('./LAVD_func')

%% set name
expName = 'JHISOdirect_32x32x16'


%% dataset definition
authkey = 'edu.jhu.pha.turbulence.testing-201406';
dataset = 'isotropic1024coarse'; % check others datasets names online


%% interpolation schemes for the JH dataset downloads
% ---- Temporal Interpolation Options ----
TInt = 'PCHIP'; % 'None', 'PCHIP'

% ---- Spatial Interpolation Flags for getVelocity & getVelocityAndPressure ----
Lag = 'Lag4'; % 'None', 'Lag4', 'Lag6', 'Lag8' 

% ---- Spatial Differentiation & Interpolation Flags for getVelocityGradient & getPressureGradient ----
FD  = 'Fd4Lag4'; % 'None_Fd4' 'None_Fd6' 'None_Fd8' 'Fd4Lag4'


%% params definition (to be setted)
% dt = 0.002; % from JH dataset readme
dt = 0.002; % from JH dataset readme
dd = 2.0*pi/1024; % from JH dataset readme
npoints_max = 128*128; % from tests on the JH readme maximum number of points that can be downloaded at a single call

% initial and final time --> to set 
t0 = 0.002; % dt*randi(1024,1);
nt = 3; % number of timestep to download
tend = t0 + (nt-1)*dt;
lagDt = dt; 

% box from which ectract the data --> to set
ncube = 64; % size of the cube to download
limx = [pi (pi+32*dd)];
limy = [pi (pi+32*dd)];
limz = [pi (pi+8*dd)];

% time vectors
tv = t0:dt:tend;
indSave = 1:1:length(tv);


% space vectors
xv  = limx(1):dd:limx(2);
yv  = limy(1):dd:limy(2);
zv  = limz(1):dd:limz(2);

nx = length(xv);
ny = length(yv);
nz = length(zv);
npoints = nx*ny*nz;

[xm ym zm]  = meshgrid(xv,yv,zv);


%% LAVD calculation from JH dataset
[VMatrixv] = getJhLAVD(xm,ym,zm,tv,indSave,npoints_max,authkey,dataset,lagDt,Lag,FD,TInt);


%% save data structure
LAVD.VMatrix = VMatrixv(:,:,:,end);
LAVD.xpi = xm; LAVD.ypi = ym; LAVD.zpi= zm; 
LAVD.textr = tv(2:end);

% save(['./data/LAVDs/LAVD_' expName '.mat'],'LAVD')


%% control plot
figure('Renderer','openGL') % isosurface of LAVD
prcv = linspace(70,95,3); % percentage of LAVD at which you want to display isosurface
clmp = colormap('lines'); % colormap 
alpha = linspace(0.1,0.8,length(prcv)); % transparency decrese with LAVD level
for i = 1:length(prcv)
    p = patch(isosurface(LAVD.xpi,LAVD.ypi,LAVD.zpi,LAVD.VMatrix,prctile(LAVD.VMatrix(:),prcv(i))));
    isonormals(LAVD.xpi,LAVD.ypi,LAVD.zpi,LAVD.VMatrix,p)
    set(p,'FaceColor',clmp(i,:));
    set(p,'EdgeColor','none');
    set(p,'FaceAlpha',alpha(i));
    hold on
end
box on; grid off
view(30,60)
daspect([1 1 1])
camproj perspective
camlight('left')


figure % a random slice of LAVD
imagesc(LAVD.VMatrix(:,:,randi(size(LAVD.VMatrix,3))))
set(gca,'xtick',[]); set(gca,'xticklabel',[])
set(gca,'ytick',[]); set(gca,'yticklabel',[])
print(gcf,fullfile(['./data/LAVDs/LAVDSlice_' expName '.png']),'-dpng','-r300');

