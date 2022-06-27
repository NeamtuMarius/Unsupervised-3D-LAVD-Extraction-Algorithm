%% main_EulerianLAVDCalculator
% This main computes the LAVD starting from flow velocity fields

%% housekeeping
clear all
close all


%% addpath
addpath('./LAVD_func/');


%% set name
expName = 'JHISO';


%% data loading
load(['./data/uFields/eulDATA_' expName '.mat'])

xq = eulDATA.xm;
yq = eulDATA.ym;
zq = eulDATA.zm;

x = unique(xq);
y = unique(yq);
z = unique(zq);

nx = length(x);
ny = length(y);
nz = length(z);

dx = mean(diff(x));
dy = mean(diff(y));
dz = mean(diff(z));

limx = [min(unique(xq)) max(unique(xq))];
limy = [min(unique(yq)) max(unique(yq))];
limz = [min(unique(zq)) max(unique(zq))];

Dx = diff(limx);
Dy = diff(limy);
Dz = diff(limz);

UT = eulDATA.u;
VT = eulDATA.v;
WT = eulDATA.w;

tspan = eulDATA.tv - eulDATA.tv(1);
dt = mean(diff(tspan));

clear eulDATA


%% evaluate vorticity on the whole volume (needed to evaluate voricity along trajectories)
[Curlx_T Curly_T Curlz_T] = omegaCalc(UT,VT,WT,dx,dy,dz,dt);


%% LAVD calculation
% definition of extraction time (must be into tspan but must not necessary belong to the vector)
t0 = 0; % initial time for extraction --
T = max(tspan); % time span (must be < of max(tspan) -- T = 0.035
dt_extr = dt/2; % time increment for extraction --
t_extr = t0:dt_extr:(t0+T); % time vector for extrtaction
limxp = mean(limx) + [-Dx/3 Dx/3];
limyp = mean(limy) + [-Dy/3 Dy/3];
limzp = mean(limz) + [-Dz/3 Dz/3];
nxp = nx;
nyp = ny;
nzp = nz;
options = odeset('RelTol',1e-3,'AbsTol',1e-3); % ODE solver options --

% def initial position for integration (must be a sub volume not to have trajs that exit the volume, must not necessary belong to the grid)
xpv = linspace(limxp(1),limxp(2),nxp);
ypv = linspace(limyp(1),limyp(2),nyp); 
zpv = linspace(limzp(1),limzp(2),nzp); 
[xpi,ypi,zpi] = meshgrid(xpv,ypv,zpv); % this is a grif, I will give the vectors and reshape afterwards

% LAVD evaluation
% trajectory integration (save also vorticity along trajectories)
[xpp_t,ypp_t,zpp_t,Curlxp_t,Curlyp_t,Curlzp_t] = IntTrajT(xpi(:),ypi(:),zpi(:),tspan,options,UT,VT,WT,x,y,z,Curlx_T,Curly_T,Curlz_T,t_extr);

% reshaping to have matrices
Curlxp_t = reshape(Curlxp_t',[size(xpi,1) size(xpi,2) size(xpi,3) size(xpp_t,1)]);
Curlyp_t = reshape(Curlyp_t',[size(xpi,1) size(xpi,2) size(xpi,3) size(xpp_t,1)]);
Curlzp_t = reshape(Curlzp_t',[size(xpi,1) size(xpi,2) size(xpi,3) size(xpp_t,1)]);
  
% space average vorticity for each instant
Curlx_avg_t = nanmean(nanmean(nanmean(Curlxp_t)));
Curly_avg_t = nanmean(nanmean(nanmean(Curlyp_t)));
Curlz_avg_t = nanmean(nanmean(nanmean(Curlzp_t)));


% LVD (not jet integrate)
LVD = sqrt((Curlxp_t - Curlx_avg_t).^2 + (Curlyp_t - Curly_avg_t).^2 + (Curlzp_t - Curlz_avg_t).^2);


% LAVD (time integral average of the LVD)
VMatrix = trapz(LVD*dt,4)/T;


%% saving LAVD matrix and cohordinate in a structure
LAVD.VMatrix = VMatrix;
LAVD.xpi = xpi; LAVD.ypi = ypi; LAVD.zpi= zpi; 
LAVD.textr = t_extr;

save(['./data/LAVDs/LAVD_' expName '.mat'],'LAVD')


%% control plot
figure('Renderer','openGL') % isosurface of LAVD
prcv = linspace(75,99,3); % percentage of LAVD at which you want to display isosurface
clmp = colormap('lines'); % colormap 
alpha = linspace(0.1,0.8,length(prcv)); % transparency decrese with LAVD level
for i = 1:length(prcv)
    p = patch(isosurface(LAVD.xpi,LAVD.ypi,LAVD.zpi,LAVD.VMatrix,prctile(LAVD.VMatrix(:),prcv(i))));
    isonormals(LAVD.xpi,LAVD.ypi,LAVD.zpi,LAVD.VMatrix,p)
    set(p,'FaceColor',clmp(i,:));
    set(p,'EdgeColor','none');
    set(p,'FaceAlpha',alpha(i));
    hold on
    box on; grid off
view(30,60)
daspect([1 1 1])
camproj perspective
camlight('left')
    pause(.1)
end



