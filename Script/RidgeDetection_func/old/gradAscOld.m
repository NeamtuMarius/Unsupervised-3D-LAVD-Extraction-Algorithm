function [xpp_t ypp_t zpp_t Vn_t] = gradAsc(xpi,ypi,zpi,stepNum,xm,ym,zm,Vm,options,cpFlag)

% normalizing data
% normalize features (preserving the aspect ratio of the domain)
x = unique(xm);
y = unique(ym);
z = unique(zm);

dx = mean(diff(x));
dy = mean(diff(y));
dz = mean(diff(z));

Dx = max(x) - min(x);
Dy = max(y) - min(y);
Dx = max(z) - min(z);

% tspan = tspan*min([dx dy dz])*100;

xn = (x - min(x))/Dx;
yn = (y - min(y))/Dx;
zn = (z - min(z))/Dx;

xpin = (xpi - min(x))/Dx;
ypin = (ypi - min(y))/Dx;
zpin = (zpi - min(z))/Dx;

% normalizzo Vm
Vn = (Vm - mean(Vm(:)))/std(Vm(:));

% valuto gradiente
[Vnx Vny Vnz] = gradient(Vn);

% valuto hessiano
[Vnxx Vnxy Vnxz] = gradient(Vnx);
[Vnyx Vnyy Vnyz] = gradient(Vnx);
[Vnzx Vnzy Vnzz] = gradient(Vnx);

% valuto integration time step & tspan
[e1,e2,e3] = EigEval3D(Vnxx,Vnxy,Vnxz,Vnyx,Vnyy,Vnyz,Vnzx,Vnzy,Vnzz);
L = max(abs([e1(:); e2(:); e3(:)]));
dt = 1/(500*L);
tspan = 0:dt:dt*stepNum;
% tspan = linspace(0,1,400)


% Defining the gradient field interpolator
Vn_interp = griddedInterpolant({xn,yn,zn,tspan},permute(repmat(Vn,1,1,1,length(tspan)),[2,1,3,4]),'cubic','none');
Vnx_interp = griddedInterpolant({xn,yn,zn,tspan},permute(repmat(Vnx,1,1,1,length(tspan)),[2,1,3,4]),'cubic','none');
Vny_interp = griddedInterpolant({xn,yn,zn,tspan},permute(repmat(Vny,1,1,1,length(tspan)),[2,1,3,4]),'cubic','none');
Vnz_interp = griddedInterpolant({xn,yn,zn,tspan},permute(repmat(Vnz,1,1,1,length(tspan)),[2,1,3,4]),'cubic','none');

u = @(xn,yn,zn,t) Vnx_interp(xn,yn,zn,t);
v = @(xn,yn,zn,t) Vny_interp(xn,yn,zn,t);
w = @(xn,yn,zn,t) Vnz_interp(xn,yn,zn,t);


% Computing gradient paths and Vn along them
[~,F] = ode45(@ODEfun,tspan,[xpin; ypin; zpin],options,u,v,w);

xppn_t = F(:,1:end/3);
yppn_t = F(:,end/3+1:2*end/3);
zppn_t = F(:,2*end/3+1:end);

for i = 1:length(tspan)
    Vn_t(i,:) = Vn_interp(xppn_t(i,:),yppn_t(i,:),zppn_t(i,:),repmat(tspan(i),1,size(xppn_t,2)));
    
end


% backscaling of the data
xpp_t = xppn_t*Dx + min(x);
ypp_t = yppn_t*Dx + min(y);
zpp_t = zppn_t*Dx + min(z);



%% control figure
if cpFlag == 1
    figure % control figure for gradient ascend
    for i = 1:length(tspan)
        plot3(xpp_t(i,:),ypp_t(i,:),zpp_t(i,:),'k.');
        xlim([min(x) max(x)]); ylim([min(y) max(y)]); zlim([min(z) max(z)]);
        title('gradient climbing control figure')
        camproj('perspective')
        grid off
        box on
        view(3)
        axis equal
        pause(.01)
        
    end
    
end


end