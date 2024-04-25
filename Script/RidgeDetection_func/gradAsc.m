function [yp xp zp] = gradAsc(xp0,yp0,zp0,stepNum,xm,ym,zm,Vm,relDiff,cpFlag)
% writing feature space
x = unique(xm); y = unique(ym); z = unique(zm);
dx = mean(diff(x)); dy = mean(diff(y)); dz = mean(diff(z));
Dx = max(x) - min(x); Dy = max(y) - min(y); Dz = max(z) - min(z);

% evaluate gradientv & hessian and its eigs
[Vmx Vmy Vmz] = gradient(Vm);
[Vmxx Vmxy Vmxz] = gradient(Vmx,dx,dy,dz);
[Vmyx Vmyy Vmyz] = gradient(Vmx,dx,dy,dz);
[Vmzx Vmzy Vmzz] = gradient(Vmx,dx,dy,dz);
[e1,e2,e3] = EigEval3D(Vmxx,Vmxy,Vmxz,Vmyx,Vmyy,Vmyz,Vmzx,Vmzy,Vmzz);

% evaluate timestep by imposing the Lipscitz condition
L = max(abs([e1(:); e2(:); e3(:)]));
dt = 1/(10*L);
tspan = 0:dt:dt*stepNum;

% Defining the gradient field interpolator
disp('building interpolation functions for Vm gradient;')
Vm_interp = griddedInterpolant({x,y,z},permute(Vm,[2,1,3]),'spline','none');
Vmx_interp = griddedInterpolant({x,y,z},permute(Vmx,[2,1,3]),'spline','none');
Vmy_interp = griddedInterpolant({x,y,z},permute(Vmy,[2,1,3]),'spline','none');
Vmz_interp = griddedInterpolant({x,y,z},permute(Vmz,[2,1,3]),'spline','none');

% ascending the gradient
if cpFlag == 1; figure; end % turning on the figure

xp = xp0;
yp = yp0;
zp = zp0;



divn = stepNum/100;
j = 1;
disp('climbing the gradient')
for i = 1:stepNum % cicle of grad asc
    % saving old position
    xp0 = xp;
    yp0 = yp;
    zp0 = zp; 
    
    % getting new position
    xp = xp + dt*Vmx_interp(xp,yp,zp);
    yp = yp + dt*Vmy_interp(xp,yp,zp);
    zp = zp + dt*Vmz_interp(xp,yp,zp);
    
    % evaluating nloss & nlossGrad 
    if ~mod(i,divn) == 1 % only every divn instants
    n0 = histcnd(xp0,yp0,zp0,x,y,z);
    n = histcnd(xp,yp,zp,x,y,z);
    dn = (n-n0)>0;
    nloss(j,1) = nansum(dn(:));
    nlossthr = max(nloss)*relDiff;
    nlossGrad = gradient(nloss);
    
    if nlossthr > nloss(j,1)
        return
    end
    
    
    j = j+1;
    
    end
    
    % control plot
    if cpFlag == 1 & ~mod(i,divn) == 1 
        tiledlayout(2,1)
        nexttile % value
        plot(nloss,'b.-','MarkerSize',10); hold on
        plot([0 stepNum/divn],[nlossthr nlossthr],'r-')
        legend('loss','tolerance')
        xlim([0 stepNum/divn])
        title('value')
        nexttile % increments
        plot(nlossGrad,'b.-','MarkerSize',10)
        xlim([0 stepNum/divn])
        title('increment')
        pause(.000001)
    end
    
    
    
end

disp('done!')


end
