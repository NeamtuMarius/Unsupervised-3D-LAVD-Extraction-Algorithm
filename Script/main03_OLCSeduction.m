%% housekeeping
clear all
close all
clc


%% addpath
addpath('./OLCSs_func');


%% setting input
expName = 'JHISO';


%% loading data
load(['./data/LAVDs/LAVD_' expName '.mat'])
load(['./data/Ridges/RidgesStreamlines_' expName '.mat'])

Vm  = LAVD.VMatrix;
xm = LAVD.xpi; ym = LAVD.ypi; zm = LAVD.zpi;
xv = unique(xm)'; yv = unique(ym)'; zv = unique(zm)';
dx = mean(diff(xv)); dy = mean(diff(yv)); dz = mean(diff(zv));
dd = [dx dy dz];

lines = ridges.lines;

%% getting the planes normal to the local tangent of the ridges
tic
for i = 1:length(lines) 
    
    line = lines{i};
    line(isnan(line(:,1)),:) = []; % remove nans 
    lines{i} = line; % remove the nans also from the lines
    for j = 1:size(line,1)

        pixel_coord = [1 1 1] + ((line(j,1:3)-[xm(1) ym(1) zm(1)])./([xm(end) ym(end) zm(end)]-[xm(1) ym(1) zm(1)])).*([size(xm,1) size(xm,2) size(xm,3)] - [1 1 1]);
        

        [VSlices{i}{j}, xSlices{i}{j}, ySlices{i}{j}, zSlices{i}{j}] = obliqueslice(permute(Vm,[2 1 3]),pixel_coord,line(j,4:6),'OutputSize','full');
        xSlices{i}{j} = xm(1) + (xm(end) - xm(1))*(xSlices{i}{j} - 1)/(size(xm,1)-1);
        ySlices{i}{j} = xm(1) + (ym(end) - ym(1))*(ySlices{i}{j} - 1)/(size(ym,1)-1);
        zSlices{i}{j} = zm(1) + (zm(end) - zm(1))*(zSlices{i}{j} - 1)/(size(zm,1)-1);

    end

end
toc


%% computation of ountermost almost convex-contours 
tic
MinLength=dd(1)*5;      % Minimum length of the contours defined as 5 times the flow-field resolution 
DeficiencyThresh=1e-2;  % Parameter for defining closed-curves convexity
Nval=50;                % Number of Isocontours to be extracted
 
for i = 1:length(lines)
    line = lines{i};
        for j = 1:5:size(line,1)

            VSlice = VSlices{i}{j};
            xSlice = xSlices{i}{j};
            ySlice = ySlices{i}{j};
            zSlice = zSlices{i}{j};
            
            bnd{i}{j} = ContExt_incl(Vm,xm,ym,zm,xSlice,ySlice,zSlice,Nval,MinLength,DeficiencyThresh,lines{i}(j,[1 2 3]));
            

        end
end
toc

%% Control figure
figure
for i = 1:length(lines)
    line = lines{i};
    
    col=[rand(1) rand(1) rand(1)];
    plot3(line(:,1),line(:,2),line(:,3),'Color',col,'LineWidth',3); hold on

    for j = 1:5:size(line,1)
        if length(bnd{i}{j})>0

            plot3(bnd{i}{j}.xc{1},bnd{i}{j}.yc{1},bnd{i}{j}.zc{1},'Color',col,'LineWidth',2)
            pause(1e-16)
        end
    end
end

