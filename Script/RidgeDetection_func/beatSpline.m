function [linef] = beatSpline(line,minlr,pl)
%beatSpline beat luthi like splining. It's quite inefficiently implemented,
%would require a bit of optimization
if size(line,1) > pl*2

t = (1:size(line,1))';
dt = 1;

plv = pl*ones(size(t)); 
plv(1:pl) = (1:pl)-1;
plv(end-pl+1:end) = fliplr((1:pl)-1);

linef = nan(size(line,1),3);
tangent = nan(size(line,1),3);

warning('off')
for i = pl:length(t)-pl
    
    pli = plv(i);
    
    ti = t(i);

    A = [ones(pli*2+1,1), (ti+(-pli:pli)'*dt), (ti+(-pli:pli)'*dt).^2, (ti+(-pli:pli)'*dt).^3];

    xi = line(i-pli:i+pli,1);
    yi = line(i-pli:i+pli,2);
    zi = line(i-pli:i+pli,3);

    cxi = (A'*xi)'*inv(A'*A);
    cyi = (A'*yi)'*inv(A'*A);
    czi = (A'*zi)'*inv(A'*A);

    xif = cxi(1) + cxi(2)*ti + cxi(3)*ti^2 + cxi(4)*ti^3;
    yif = cyi(1) + cyi(2)*ti + cyi(3)*ti^2 + cyi(4)*ti^3;
    zif = czi(1) + czi(2)*ti + czi(3)*ti^2 + czi(4)*ti^3;

    txif = cxi(2) + 2*cxi(3)*ti + 3*cxi(4)*ti^2;
    tyif = cyi(2) + 2*cyi(3)*ti + 3*cyi(4)*ti^2;
    tzif = czi(2) + 2*czi(3)*ti + 3*czi(4)*ti^2;
    
    linef(i,1) = xif;
    linef(i,2) = yif;
    linef(i,3) = zif;

    tangent(i,1) = txif./sqrt(txif.^2 + tyif.^2 + tzif.^2);
    tangent(i,2) = tyif./sqrt(txif.^2 + tyif.^2 + tzif.^2);
    tangent(i,3) = tzif./sqrt(txif.^2 + tyif.^2 + tzif.^2);

end
warning('on')

else

linef = nan(size(line,1),3);
tangent = nan(size(line,1),3);

end

linef = [linef tangent];

end