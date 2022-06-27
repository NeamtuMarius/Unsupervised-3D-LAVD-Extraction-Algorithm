function [VMatrix] = getJhLAVD(xm,ym,zm,tv,indSave,npoints_max,authkey,dataset,lagDt,Lag,FD,TInt)

nx = size(xm,1);
ny = size(xm,2);
nz = size(xm,3);
npoints = nx*ny*nz;
nt = length(tv);
t0 = tv(1);
dt = mean(diff(tv));
tvSave = tv(indSave);

points = [xm(:)'; ym(:)'; zm(:)'];
[intv] = n2intv(npoints,npoints_max);
npointsi = diff(intv,1) + 1;
LVDint = zeros(nx,ny,nz);
T = tv(end) - tv(2);
LVDint = zeros(nx,ny,nz);
jt = 1;

for it = 1:nt
    tic;
    for i = 1:size(intv,2)
        if it == 1
            pos(:,intv(1,i):intv(2,i)) = points(:,intv(1,i):intv(2,i)); % position
            duAll(:,intv(1,i):intv(2,i)) = getVelocityGradient(authkey,dataset,tv(it),FD,TInt,npointsi(i),pos(:,intv(1,i):intv(2,i)));
        else
            pos(:,intv(1,i):intv(2,i)) = getPosition(authkey, dataset, tv(it-1), tv(it), lagDt, Lag, npointsi(i),pos(:,intv(1,i):intv(2,i))); % position
            duAll(:,intv(1,i):intv(2,i)) = getVelocityGradient(authkey,dataset,tv(it),FD,TInt,npointsi(i),pos(:,intv(1,i):intv(2,i)));
        end

        fprintf('\nExtracted %.2f per cent of the whole time and %.2f per cent of the whole volume\n',tv(it)/tv(end)*100,i/size(intv,2)*100)
        
    end
    
    % shaping data in 4D matrices
    xpos = reshape(pos(1,:),nx,ny,nz);
    ypos = reshape(pos(2,:),nx,ny,nz);
    zpos = reshape(pos(3,:),nx,ny,nz);
    ux = reshape(duAll(1,:),nx,ny,nz); uy = reshape(duAll(2,:),nx,ny,nz); uz = reshape(duAll(3,:),nx,ny,nz);
    vx = reshape(duAll(4,:),nx,ny,nz); vy = reshape(duAll(5,:),nx,ny,nz); vz = reshape(duAll(6,:),nx,ny,nz);
    wx = reshape(duAll(7,:),nx,ny,nz); wy = reshape(duAll(8,:),nx,ny,nz); wz = reshape(duAll(9,:),nx,ny,nz);
    
    % vorticity
    Curlxp_t = wy - vz;
    Curlyp_t = uz - wx;
    Curlzp_t = vx - uy;
    
    % spatial average of vorticity
    Curlx_avg_t = nanmean(nanmean(nanmean(Curlxp_t)));
    Curly_avg_t = nanmean(nanmean(nanmean(Curlyp_t)));
    Curlz_avg_t = nanmean(nanmean(nanmean(Curlzp_t)));
    
    % Lag vorticity dev
    LVD = sqrt((Curlxp_t - Curlx_avg_t).^2 + (Curlyp_t - Curly_avg_t).^2 + (Curlzp_t - Curlz_avg_t).^2);
    
    LVDint = LVDint + dt*LVD;
    
    if any(it == indSave)
        TT = tvSave(jt) - tvSave(1);
        VMatrix(:,:,:,jt) = LVDint/TT;
        jt = jt+1;
        
    end
    
    timeval(it) = toc;
    % fprintf('\nExtracted %.2f per cent of the whole time in %.2f sec\n',tv(it)/tv(end)*100,timeval(it))
    % fprintf('\nTotal elapsed time: %.2f sec of (expected) %.2f sec\n',cumsum(timeval),mean(timeval(it))*nt)
end

end