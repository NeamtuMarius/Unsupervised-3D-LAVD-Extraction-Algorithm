%% References:
% [1]: G. Haller, A. Hadjighasem, M. Farazamand & F. Huhn,Vm
%      Defining coherent vortces objectively from the vorticity. submitted (2015)
% [2]: G. Haller, Dynamically consistent rotation and stretch tensors for
%      finite continuum deformation. submitted (2015).

function bnd = ContExt_incl(Vm,xm,ym,zm,xSlice,ySlice,zSlice,Nval,MinLength,DeficiencyThresh,P)

% Curves Extraction in inclined planes

Vm=permute(Vm,[2 1 3]);
Nct=linspace(prctile(Vm(:),50),prctile(Vm(:),99),Nval);


        figure('visible','off')
        contourslice(xm,ym,zm,Vm,xSlice,ySlice,zSlice,Nct); 
%         [M,c] = contour3(xm,ym,zm,Vm,xSlice,ySlice,zSlice,Nct);
        h=gcf;
        axesObjs = get(h, 'Children');  %axes handles
        dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
        xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
        ydata = get(dataObjs, 'YData');
        zdata = get(dataObjs, 'ZData');
        cdata=get(dataObjs,'CData');
        
% defining the rotation matrix        
        xcum=[];ycum=[];zcum=[];
        for i=1:length(xdata)
            xcum=[xcum xdata{i}(1:end-1)'];
            ycum=[ycum ydata{i}(1:end-1)'];
            zcum=[zcum zdata{i}(1:end-1)'];
        end
        
        xint0 = xcum-nanmean(xcum) ; 
        yint0 = ycum-nanmean(ycum) ;
        zint0 = zcum-nanmean(zcum) ;
        
        PCA=pca([xint0' yint0' zint0']);
        e1=PCA(:,1)'; e2=PCA(:,2)' ;e3=PCA(:,3)';  % 3 principal vector(3 eigenvector) of "input data"
        % n1=[1 0 0]  ; n2=[0 1 0]   ; n3=[0 0 1];   % 3 unit vector Ox,Oy,Oz
        %  transformation matrix from "e" space to "n" space
        R=[e1;e2;e3];  % rotation matrix , match with xyz (e1//n1, e2//n2, e3//n3) 
        % Finding the new rotate data
%         newdata=(R*[xint0' yint0' zint0']')';%new data corresponding to P1 coordinat
        
        
for j=1:length(xdata)
    app=(R*[xdata{j}(1:end-1)-nanmean(xcum) ydata{j}(1:end-1)-nanmean(ycum) zdata{j}(1:end-1)-nanmean(zcum)]')';
    xdata{j}=app(:,1);
    ydata{j}=app(:,2);
    zdata{j}=app(:,3);
end

Pnew=(R*[P(1)-nanmean(xcum) P(2)-nanmean(ycum) P(3)-nanmean(zcum)]')';


sprintf('... %3d contours are extracted.',numel(xdata))

bnd = struct('xc',[],'yc',[],'zc',[],'cval',[],'xp',[],'yp',[],'zp',[],'valp',[]);
for k=1:numel(xdata)
    % check if the contour is closed
    if ( xdata{k}(1)==xdata{k}(end) ) && ( ydata{k}(1)==ydata{k}(end) )
        Length = sum( sqrt(diff(xdata{k}).^2+diff(ydata{k}).^2) );
        if Length>MinLength
            % Check if the contour is convex
            if IsContourConvex( xdata{k}, ydata{k},DeficiencyThresh) 
                bnd.xc{end+1} = xdata{k};
                bnd.yc{end+1} = ydata{k};
                bnd.zc{end+1} = zdata{k}; 
                bnd.cval(end+1,1) = cdata{k}(1,1);
            end
        end
    end
end
Nct_filt1 = numel(bnd.xc);
sprintf('... %3d contours are closed and convex.',Nct_filt1)
%% Step2: Finding local maxima of the LAVD field.

if Nct_filt1>0

%- Keeping closed contours that encompass only "1" local maximum.
InMax = cellfun(@(x,y) inpolygon(Pnew(1),Pnew(2),x,y),bnd.xc,bnd.yc,'UniformoutPut',false);
InMax = cat(1,InMax{:});       % rows --> closed contours & columns --> local maxima
N_InMax = sum(InMax,2);        % Number of local maxima in each contour
indEliminate_1 = N_InMax~=1;
bnd.xc(indEliminate_1)   = [];
bnd.yc(indEliminate_1)   = [];
bnd.zc(indEliminate_1)   = [];
bnd.cval(indEliminate_1) = [];


%% Step3: Selecting the outermost contour for each nested family of contours
%- selecting the first point of each contour as a query point
Nct_filt2 = numel(bnd.xc);
if Nct_filt2~=0
    xq = cellfun(@(x) x(1),bnd.xc,'UniformoutPut',true);  
    yq = cellfun(@(x) x(1),bnd.yc,'UniformoutPut',true);  
end
%%
indEliminate_2 = false(Nct_filt2,1);
for ii=1:Nct_filt2
    in = inpolygon(xq,yq,bnd.xc{ii},bnd.yc{ii});  in(ii)=0;
    indEliminate_2(in) = 1;
end
bnd.xc(indEliminate_2)   = [];
bnd.yc(indEliminate_2)   = [];
bnd.zc(indEliminate_2)   = [];
bnd.cval(indEliminate_2) = [];

if  numel(bnd.xc)==1
    
inv_data=(inv(R)*[bnd.xc{1} bnd.yc{1} bnd.zc{1}]')';
inv_data(:,1)=inv_data(:,1)+nanmean(xcum);
inv_data(:,2)=inv_data(:,2)+nanmean(ycum);
inv_data(:,3)=inv_data(:,3)+nanmean(zcum);

bnd.xc{1}=inv_data(:,1);
bnd.yc{1}=inv_data(:,2);
bnd.zc{1}=inv_data(:,3);
else
  bnd=[];  
end

else
bnd=[];
end

clear inv_data xdata ydata zdata xcum ycum zcum xint0 yint0 zint0 app cdata InMax...
    xm ym zm Vm xSlice ySlice zSlice axesObjs dataObjs
close all
end
