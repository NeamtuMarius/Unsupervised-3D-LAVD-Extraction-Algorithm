function [filteredLine] = line2filteredLine(line,lpf,minlr)
%line2filteredLine apply a low pass filter and discarge the first and last
%minlr points
%   

x = line(:,1);
y = line(:,2);
z = line(:,3);

dx = [diff(x); diff(x(end-1:end))];
dy = [diff(y); diff(y(end-1:end))];
dz = [diff(z); diff(z(end-1:end))];

ds = sqrt(dx.^2 + dy.^2 + dz.^2);
s = cumsum(ds);
dsnom = nanmean(ds);
fsnom = 1/dsnom;

[xrs, srs] = resample(x,s,fsnom,'linear');
[yrs, srs] = resample(y,s,fsnom,'linear');
[zrs, srs] = resample(z,s,fsnom,'linear');

xf = lowpass(xrs,lpf); xf([1:minlr end-minlr+1:end]) = []; 
yf = lowpass(yrs,lpf); yf([1:minlr end-minlr+1:end]) = [];
zf = lowpass(zrs,lpf); zf([1:minlr end-minlr+1:end]) = [];

filteredLine = [xf yf zf];


end

