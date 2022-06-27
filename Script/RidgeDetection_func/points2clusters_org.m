function [clustersLong] = points2clusters_org(points,distmin,prc)
%points2clusters transform the point clouds into clusters of point  belonging  to lines
% of minimal length and minimal distance between points
%   distmin: minimum euclidean distance between 2 points to belong to the same line
%   lenmin: minimum length of a line two be so (in terms of vector length, TOO BE MODIFIED)

pointsColor = double(pcsegdist(pointCloud(points),distmin));
labels = unique(pointsColor);

for i = 1:length(labels)
    cluster = points(pointsColor == labels(i),:);
    clusters{i} = cluster;
    clusterLength(i,1) = size(clusters{i},1);
    
end

minLength = prctile(clusterLength,prc);
indLong = find(clusterLength >= minLength);

clustersLong = clusters(indLong);


end

