function lines = clusters2lines(clusters,dsmin)
%clusters2lines Summary of this function goes here
%   Detailed explanation goes here

disp('building the lines...     ')

clusterColors = unique(clusters(:,end));

for i = 1:length(clusterColors)
    
    line = cluster2line(clusters(clusters(:,end) == clusterColors(i),:));
    
    xl = line(:,1);
    yl = line(:,2);
    zl = line(:,2);

    dxl = gradient(xl);
    dyl = gradient(yl);
    dzl = gradient(zl);
    
    dsl = sqrt(dxl.^2 + dyl.^2 + dzl.^2);
    indbad = find(dsl>dsmin);

    dsl_new = dsl; dsl_new(indbad) = nan;
    
    
    lines{i} = line; lines{i}(indbad,:) = nan;
    fprintf('\b\b\b\b%3.0f%%', (100*(i/(length(clusterColors)))));
    
end


end

