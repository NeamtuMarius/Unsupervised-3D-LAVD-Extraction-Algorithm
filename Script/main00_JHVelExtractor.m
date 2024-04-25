%% main_JHTurbExtractor
% extract flow velocity data data from jht database without any size limitation
% (need to download "turbmat-master" from "https://github.com/idies/turbmat" and add the path of JHT database matlab functions)

%% housekeeping
clear all;
close all;

%% addpath of the JH turb dataset matlab func folder
addpath('./turbmat-master'); % JH functions
addpath('./LAVD_func/'); 

%% dataset definition
authkey = 'edu.jhu.pha.turbulence.testing-201406';
dataset = 'isotropic1024coarse'; % check others datasets names online


%% interpolation schemes for the JH dataset downloads
% ---- Temporal Interpolation Options ----
TInt = 'PCHIP'; % 'None', 'PCHIP'

% ---- Spatial Interpolation Flags for getVelocity & getVelocityAndPressure ----
Lag = 'Lag4'; % 'None', 'Lag4', 'Lag6', 'Lag8' 


%% params definition (to be setted)
dt = 0.002; % from JH dataset readme
dd = 2.0*pi/1024; % from JH dataset readme
npoints_max = 64*64; % from tests on the JH readme maximum number of points that can be downloaded at a single call

% initial and final time --> to set 
t0 = dt; % dt*randi(1024,1);
nt = 20; % number of timestep to download
tend = t0 + (nt-1)*dt;

% box from which ectract the data --> to set
ncube = 128; % size of the cube to download
limx = [pi (pi+(ncube-1)*dd)];
limy = limx;
limz = limx;

% time vectors
tv = t0:dt:tend;

% space vectors
xv  = limx(1):dd:limx(2);
yv  = limy(1):dd:limy(2);
zv  = limz(1):dd:limz(2);

nx = length(xv);
ny = length(yv);
nz = length(zv);
npoints = nx*ny*nz;

[xm ym zm]  = meshgrid(xv,yv,zv);


%% eul data extraction extradction
points = [xm(:)'; ym(:)'; zm(:)'];
[intv] = n2intv(npoints,npoints_max);
npointsi = diff(intv,1) + 1;
u = zeros(nx,ny,nz,nt);
v = zeros(nx,ny,nz,nt); 
w = zeros(nx,ny,nz,nt);
for it = 1:nt
    
    for i = 1:size(intv,2)
        vel(:,intv(1,i):intv(2,i)) = getVelocity(authkey,dataset,tv(it),Lag,TInt,npointsi(i),points(:,intv(1,i):intv(2,i)));
        
        fprintf('\nTime %.2f of %.2f, extracted %.2f per cent of the whole volume\n',tv(it),tv(end),i/size(intv,2)*100)
        
    end
    
    % shaping data in 4D matrices
    u(:,:,:,it) = reshape(vel(1,:),nx,ny,nz);
    v(:,:,:,it) = reshape(vel(2,:),nx,ny,nz);
    w(:,:,:,it) = reshape(vel(3,:),nx,ny,nz);
    

end


%% dataset saving
% create a structure
eulDATA.xm = xm; eulDATA.ym = ym; eulDATA.zm = zm; eulDATA.tv = tv; % grids
eulDATA.u = u; eulDATA.v = v; eulDATA.w = w; % velocity

% saving
save('./data/uFields/eulDATA_JHISO.mat','eulDATA') % save manually if it does not work


%% control plot
figure
i  = 32;
for it = 1:length(eulDATA.tv)
    
    imagesc('xdata',unique(eulDATA.xm),'ydata',unique(eulDATA.ym),'cdata',sqrt(eulDATA.u(:,:,i,it).^2 + eulDATA.v(:,:,i,it).^2 + eulDATA.w(:,:,i,it).^2)); hold on
    xlim([min(eulDATA.xm(:)) max(eulDATA.xm(:))]); ylim([min(eulDATA.ym(:)) max(eulDATA.ym(:))]);
    set(gca,'ydir','normal')
    set(gca,'xtick',[]); set(gca,'xticklabel',[])
    set(gca,'ytick',[]); set(gca,'yticklabel',[])
    colormap(hot)
    pause(.001)
    
    hold off
    
end





