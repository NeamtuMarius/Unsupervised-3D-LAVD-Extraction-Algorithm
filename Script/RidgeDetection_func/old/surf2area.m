function [area] = surf2area(verts,faces)
%surf2area gets verts and faces and evaluate the area
%   Detailed explanation goes here
a = verts(faces(:, 2), :) - verts(faces(:, 1), :);
b = verts(faces(:, 3), :) - verts(faces(:, 1), :);
c = cross(a, b, 2);
area = 1/2 * sum(sqrt(sum(c.^2, 2)));
end

