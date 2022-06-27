function [sorted_line] = lineSorting(line)
%lineSorting Summary of this function goes here
%   Detailed explanation goes here
dist = pdist2(line,line);

N = size(line,1);
result = NaN(1,N);
result(1) = 1; % first point is first row in data matrix

for ii=2:N
    dist(:,result(ii-1)) = Inf;
    [~, closest_idx] = min(dist(result(ii-1),:));
    result(ii) = closest_idx;
    
end

sorted_line = line(result(1:end-1),:);


end

