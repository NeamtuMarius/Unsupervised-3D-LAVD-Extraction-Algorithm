function [xpp_t,ypp_t,zpp_t,Curlpx_t,Curlpy_t,Curlpz_t] = IntTrajT(xpi,ypi,zpi,tspan,options,UT,VT,WT,x,y,z,Curlx_T,Curly_T,Curlz_T,timeint)
tic
% Np = numel(xi);
% Nt = numel(tspan);
%% Defining the velocity field
u_interp = griddedInterpolant({x,y,z,tspan},permute(UT,[2,1,3,4]),'linear','none');
v_interp = griddedInterpolant({x,y,z,tspan},permute(VT,[2,1,3,4]),'linear','none');
w_interp = griddedInterpolant({x,y,z,tspan},permute(WT,[2,1,3,4]),'linear','none');

u = @(x,y,z,t) u_interp(x,y,z,t);
v = @(x,y,z,t) v_interp(x,y,z,t);
w = @(x,y,z,t) w_interp(x,y,z,t);
        
%% Computing Lagrangian trajectories:
[~,F] = ode45(@ODEfun,timeint,[xpi;ypi;zpi],options,u,v,w);

xpp_t = F(:,1:end/3);
ypp_t = F(:,end/3+1:2*end/3);
zpp_t = F(:,2*end/3+1:end);


%% Computing vorticity along trajectories:
Culx_interp = griddedInterpolant({x,y,z,tspan},permute(Curlx_T,[2,1,3,4]),'linear','none');
Culy_interp = griddedInterpolant({x,y,z,tspan},permute(Curly_T,[2,1,3,4]),'linear','none');
Culz_interp = griddedInterpolant({x,y,z,tspan},permute(Curlz_T,[2,1,3,4]),'linear','none');

for i = 1:length(timeint)
    Curlpx_t(i,:) = Culx_interp(xpp_t(i,:),ypp_t(i,:),zpp_t(i,:),repmat(timeint(i),1,size(xpp_t,2)));
    Curlpy_t(i,:) = Culy_interp(xpp_t(i,:),ypp_t(i,:),zpp_t(i,:),repmat(timeint(i),1,size(xpp_t,2)));
    Curlpz_t(i,:) = Culz_interp(xpp_t(i,:),ypp_t(i,:),zpp_t(i,:),repmat(timeint(i),1,size(xpp_t,2)));
    
end

end