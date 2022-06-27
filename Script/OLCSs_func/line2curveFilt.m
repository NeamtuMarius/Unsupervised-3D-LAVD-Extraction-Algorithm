function [curve] = line2curveFilt(line,ff)
%line2linetFilt get line, resample the fata on uniform curvilinear
%coordinates, evaluate t, and filter it
%   
xs = fillmissing(line(:,1),'linear');
ys = fillmissing(line(:,2),'linear');
zs = fillmissing(line(:,3),'linear');

dxs = [diff(xs); diff(xs(end-1:end))];
dys = [diff(ys); diff(ys(end-1:end))];
dzs = [diff(zs); diff(zs(end-1:end))];

ds = sqrt(dxs.^2 + dys.^2 + dzs.^2);
s = cumsum(ds);

fs = 1/mean(ds);
[xsrs, srs] = resample(xs,s,fs);
[ysrs, srs] = resample(ys,s,fs);
[zsrs, srs] = resample(zs,s,fs);

dxsrs = gradient(xsrs);
dysrs = gradient(ysrs);
dzsrs = gradient(zsrs);

dsrs = sqrt(dxsrs.^2 + dysrs.^2 + dzsrs.^2);
srs = cumsum(dsrs);
 
tx = dxsrs./dsrs;
ty = dysrs./dsrs;
tz = dzsrs./dsrs;

% ff = 0.1; % fraction for filtering

txf = medfilt1(tx,floor(length(tx)*ff));
tyf = medfilt1(ty,floor(length(tx)*ff));
tzf = medfilt1(tz,floor(length(tx)*ff));

pos = [xsrs ysrs zsrs];
s = srs;


curve.pos = [xsrs ysrs zsrs];
curve.s   = srs;
curve.t   = [txf tyf tzf];

end

