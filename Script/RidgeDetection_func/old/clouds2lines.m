function [lines] = clouds2lines(points,distmin,lenmin,cpFlag)
%points2clusters transform the point clouds into clusters of point  belonging  to lines
% of minimal length and minimal distance between points
%   distmin: minimum euclidean distance between 2 points to belong to the same line
%   lenmin: minimum length of a line two be so (in terms of vector length, TOO BE MODIFIED)


% clustering the cluouds
labels = pcsegdist(pointCloud(points),distmin) + 1;

% keeping the long ones
k = 1;
for i = 1:max(labels)
    
    line = points(labels == i,:);
    
    if size(line,1)>=lenmin
        lines{k} = ridgesOutliersDetector(line);
        k = k+1;
    end
end


if cpFlag == 1
    figure % control figure
    for i = 1:length(lines)
        plot3(lines{i}(:,1),lines{i}(:,2),lines{i}(:,3),'-','LineWidth',2); hold on
        
    end

    camproj('perspective')
    grid off
    box on
    axis equal
    view(3)
    
end

end

