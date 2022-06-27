function [sorted_line] = ridgesOutliersDetector(line)
%ridgesOutliersDetector find point not belonging to the minimum spanning
%tree longest branch
%   

stopControl = 0;
while stopControl ~= 2
    
    D = pdist2(line,line);
    g = graph(D);
    T = minspantree(g);
    
    T1 = T;
    
    idx = find(T1.degree==1);
    line1 = line; line1(idx,:) = [];
    
    
    
    stopControl = size(line,1) - size(line1,1);
    line = line1;
    
end

[sorted_line] = lineSorting(line1);


end

